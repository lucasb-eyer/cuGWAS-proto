#include "fgls.h"
#include "io.h"

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>

#include <sys/time.h>
#include <time.h>

#include "vt_user.h"

#define NUM_BUFFERS_PER_THREAD 2


typedef struct {
	FILE *X_fp;
	FILE *Y_fp;
	FILE *B_fp;

	double *X[NUM_BUFFERS_PER_THREAD];
	double *Y[NUM_BUFFERS_PER_THREAD];
	double *B[NUM_BUFFERS_PER_THREAD];

	double *h;
	double *W;
	double *Winv;
	double *XtZWinv;
	double *xtSx;

	sem_t *sem_io;
	sem_t *sem_comp_x;
	sem_t *sem_comp_y;
	
	int id;
} ooc_loops_t;


void all_sem_post(sem_t *sem, int ntimes, char *var) {
	int i;
	int value;
	for ( i = 0; i < ntimes; i++ )
	{
		/*printf("IO THREAD Posts SEM_COMP_%s[%d]\n", var, i);*/
		/*sem_getvalue(&sem[i], &value);*/
		/*printf("IO THREAD Pre  value SEM_COMP_%s[%d]: %d\n", var, i, value);*/
		sem_post( &sem[i] );
		/*sem_getvalue(&sem[i], &value);*/
		/*printf("IO THREAD Post value SEM_COMP_%s[%d]: %d\n", var, i, value);*/
	}
}
void all_sem_wait(sem_t *sem, int ntimes) {
	int i;
	int value;
	for ( i = 0; i < ntimes; i++ )
	{
		/*printf("IO THREAD Waits\n");*/
		/*sem_getvalue(&sem[i], &value);*/
		/*printf("IO THREAD Pre  value SEM_IO[%d]: %d\n", i, value);*/
 		sem_wait( &sem[i] );
		/*sem_getvalue(&sem[i], &value);*/
		/*printf("IO THREAD Post value SEM_IO[%d]: %d\n", i, value);*/
	}
}
void write_all_b( double *buff, FILE *fp, int i, int j, FGLS_eigen_t *cf) {
	int it;
	/*printf("Write all: %d %d\n", i, j);*/
	for ( it = 0; it < MIN( cf->NUM_COMPUTE_THREADS, cf->t - j ); it++ )
	{
		/*printf("IO THREAD Writes from B[%d] to File[%d]: %d\n", it * cf->x_b * cf->p, (j + it) * (cf->m * cf->p) + i * cf->p, buff[ it * cf->x_b * cf->p ]); */
        write( &buff[ it * cf->x_b * cf->p ], fp, 
				MIN( cf->x_b * cf->p, (cf->m - i) * cf->p),
				(j + it) * (cf->m * cf->p) + i * cf->p);
	}
}

void* io_thread_func(void *in) 
{
  DEF_TIMING();

  ooc_loops_t *loops_t = ( ooc_loops_t* ) in;
  FGLS_eigen_t *cf = &FGLS_eigen_config;

  double *x_cur  = loops_t->X[0];
  double *x_next = loops_t->X[1];
  double *y_cur  = loops_t->Y[0];
  double *y_next = loops_t->Y[1];
  double *b_prev = loops_t->B[0];
  double *b_cur  = loops_t->B[1];

  int m = cf->m,
	  n = cf->n,
	  p = cf->p,
	  t = cf->t,
	  x_b = cf->x_b;
  int i, j;

  BEGIN_TIMING();
  read( x_cur, loops_t->X_fp, MIN( x_b * n * p, m * p * n ), 0 );
  read( y_cur, loops_t->Y_fp, MIN( cf->NUM_COMPUTE_THREADS * n, t * n ), 0 );
  END_TIMING(cf->time->io_time);

  all_sem_post( loops_t->sem_comp_x, cf->NUM_COMPUTE_THREADS, "X" );
  all_sem_post( loops_t->sem_comp_y, cf->NUM_COMPUTE_THREADS, "Y" );
  /*printf("Both threads can proceed with Y[0]\n");*/
  for ( j = 0; j < t; j += cf->NUM_COMPUTE_THREADS ) 
  {
    swap_buffers( &y_cur, &y_next );
    
    VT_USER_START("IO_READ_Y");
    BEGIN_TIMING();
    read( y_cur, loops_t->Y_fp, 
			j + cf->NUM_COMPUTE_THREADS > t ? 0 : MIN( cf->NUM_COMPUTE_THREADS * n, (t - (j + cf->NUM_COMPUTE_THREADS)) * n ), 
			j + cf->NUM_COMPUTE_THREADS > t ? 0 : (j + cf->NUM_COMPUTE_THREADS) * n );
    END_TIMING(cf->time->io_time);
    all_sem_post( loops_t->sem_comp_y, cf->NUM_COMPUTE_THREADS, "Y" );
    VT_USER_END("IO_READ_Y");

    for ( i = 0; i < m; i += x_b ) 
	{
      swap_buffers( &x_cur, &x_next );

      VT_USER_START("IO_READ_X");
      BEGIN_TIMING();
      read( x_cur, loops_t->X_fp, 
			  (i + x_b) >= m ? MIN( x_b * n * p, m * p * n ) : MIN( x_b * n * p, (m - (i + x_b)) * p * n ), 
			  (i + x_b) >= m ? 0 : (i + x_b) * p * n );
      END_TIMING(cf->time->io_time);
      all_sem_post( loops_t->sem_comp_x, MIN( cf->NUM_COMPUTE_THREADS, (t - j) ), "X" );
	  /*printf("Both threads can proceed with Y[%d]\n", j+cf->NUM_COMPUTE_THREADS);*/
      VT_USER_END("IO_READ_X");

      BEGIN_TIMING();
      all_sem_wait( &loops_t->sem_io[0], MIN( cf->NUM_COMPUTE_THREADS, (t - j) ) );
      END_TIMING(cf->time->io_mutex_wait_time);

      swap_buffers( &b_prev, &b_cur );

      VT_USER_START("IO_WRITE_B");
      BEGIN_TIMING();
      write_all_b( b_prev, loops_t->B_fp, i, j, cf );
      END_TIMING(cf->time->io_time);
      VT_USER_END("IO_WRITE_B");
	  /*printf("B[0]: %f\n", *b_prev );*/
    }
  }
  pthread_exit(NULL);
}

void* compute_thread_func(void* in) {
  DEF_TIMING();

  ooc_loops_t *loops_t = ( ooc_loops_t* ) in;
  FGLS_eigen_t *cf = &FGLS_eigen_config;

  int m = cf->m,
	  n = cf->n,
	  p = cf->p,
	  t = cf->t,
	  x_b = cf->x_b,
	  id = loops_t->id;

  double *x_cur  =  loops_t->X[0];
  double *x_next =  loops_t->X[1];
  double *y_cur  = &loops_t->Y[0][id * n];
  double *y_next = &loops_t->Y[1][id * n];
  double *b_prev = &loops_t->B[0][id * p * x_b];
  double *b_cur  = &loops_t->B[1][id * p * x_b];
  double *h      =  loops_t->h;
  double *W      =  loops_t->W;

  double *Winv    = &loops_t->Winv[id * n];
  double *XtZWinv = &loops_t->XtZWinv[id * x_b * p * n];
  double *xtSx    = &loops_t->xtSx[id * p * p];

  double ONE = 1.0;
  double ZERO = 0.0;
  int iONE = 1;
  int iZERO = 0;
  double alpha, beta;

  int ib, i, j, k;
  int x_inc;
  int mpb;
  int mpb_real;
  int info;

  int value;
  struct timespec t0, t1;

  for (j = id; j < t; j += cf->NUM_COMPUTE_THREADS) 
  {
	  /*clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t0);*/
	  /*printf("Id: %d - self: %d\n", id, pthread_self());*/
  	alpha = h[j] * h[j];
	beta  = 1 - alpha;

	/*sem_getvalue(&loops_t->sem_comp_y[id], &value);*/
	/*printf("SEM_COMP_Y[%d]: value: %d\n", id, value);*/
    VT_USER_START("COMP_WAIT_Y");
    BEGIN_TIMING();
    sem_wait( &loops_t->sem_comp_y[id] );
    END_TIMING(cf->time->comp_mutex_wait_time);
    VT_USER_END("COMP_WAIT_Y");
	/*sem_getvalue(&loops_t->sem_comp_y[id], &value);*/
	/*printf("SEM_COMP_Y[%d]: value: %d\n", id, value);*/

    VT_USER_START("COMP_LOOP_Y_CODE");
    BEGIN_TIMING();
    /* 2) W = sqrt(alpha W - beta I)^-1 */
    // Best order? sqrt - inv
    for (k = 0; k < n; k++)
      Winv[k] = sqrt(1.0 / (alpha * W[k] + beta));
    
    /* sqrt(Winv) * ZtY */
    for (k = 0; k < n; k++)
      /*dscal_(&t, &Winv[l], &ZtY[l], &n);*/
      y_cur[k] *= Winv[k];
    END_TIMING(cf->time->compute_time);
    VT_USER_END("COMP_LOOP_Y_CODE");

    for (ib = 0; ib < m; ib += x_b) 
	{
		/*sem_getvalue(&loops_t->sem_comp_x[id], &value);*/
		/*printf("SEM_COMP_X[%d]: value: %d\n", id, value);*/
		BEGIN_TIMING();
		VT_USER_START("COMP_WAIT_X");
      sem_wait( &loops_t->sem_comp_x[id] );
		VT_USER_END("COMP_WAIT_X");
	  END_TIMING2(cf->time->comp_mutex_wait_time);
	  /*sem_getvalue(&loops_t->sem_comp_x[id], &value);*/
	  /*printf("SEM_COMP_X[%d]: value: %d\n", id, value);*/

      VT_USER_START("COMP_LOOP_X_CODE");
      BEGIN_TIMING();
      x_inc = MIN(x_b, m - ib);
	  mpb = p * x_b;
	  mpb_real = p * x_inc;

	  /*printf("X[0]: %f\n", x_cur[0]);*/
	  /*printf("Y[0]: %f\n", y_cur[0]);*/
	  /*printf("D[0]: %f\n", Winv[0]);*/
      /* X' * Z  * sqrt(Winv) */
      for (k = 0; k < n * mpb; k++) XtZWinv[k] = 0.0;
      for (k = 0; k < n; k++)
        daxpy_(&mpb_real, &Winv[k], &x_cur[k], &n, &XtZWinv[k*mpb], &iONE);
	  /*printf("W[0]: %f\n", XtZWinv[0]);*/
      
      /* 7) y = XtZWinv * y */
      dgemv_("N", &mpb_real, &n, &ONE, XtZWinv, &mpb, y_cur, &iONE, &ZERO, b_cur, &iONE);
	  /*printf("RHS[0] - id: %d - j: %d: %f\n", id, j, b_cur[0]);*/
      
	  /*printf("x_inc: %d\n", x_inc);*/
      for (i = 0; i < x_inc; i++)
      {
        /* 5) W = XtZWinv * K^T */
        dsyrk_("L", "N", &p, &n, &ONE, &XtZWinv[i*p], &mpb, &ZERO, xtSx, &p);

        /* 8) W^-1 * y */
        dposv_("L", &p, &iONE, xtSx, &p, &b_cur[i*p], &p, &info);
        if (info != 0)
        {
          fprintf(stderr, "Error executing dposv (s: %d, r: %d, c: %d, th: %d): %d\n", j, ib, i, id, info);
          exit(-1);
        }
      }
      END_TIMING(cf->time->compute_time);

      swap_buffers( &x_cur, &x_next );
      swap_buffers( &b_prev, &b_cur );
	  /*printf("COMP B[0]: %f\n", *b_prev );*/

	  /*sem_getvalue(&loops_t->sem_io[id], &value);*/
	  /*printf("COMP THREAD %d Posts SEM_IO    : value: %d\n", id, value);*/
      sem_post( &loops_t->sem_io[id] );
	  /*sem_getvalue(&loops_t->sem_io[id], &value);*/
	  /*printf("COMP THREAD %d After SEM_IO    : value: %d\n", id, value);*/
      VT_USER_END("COMP_LOOP_X_CODE");
    }
    swap_buffers( &y_cur, &y_next );
	/*gettimeofday(&t1, NULL);*/
	/*clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t1);*/
	/*printf("Iteration %3d: %10ld ms\n", j, get_diff_ns(&t0, &t1)/1000000);*/
  }
  pthread_exit(NULL);
}

int fgls_eigen( FGLS_eigen_t *cf ) {
  int rc, i;
  pthread_t io_thread;
  pthread_t *compute_threads;

  double *Phi,
		 *Z,
		 *W;
  FILE *Phi_fp,
	   *h_fp;
  ooc_loops_t loops_t;
  ooc_loops_t *loops_t_comp;

  /*sprintf(str_buf, "%s/X.in", dir);*/
  /*x_file = fopen(str_buf, "rb");*/
  /*if(!x_file) {*/
  /*printf("error opening x_file(%s)! exiting...\n", str_buf);*/
  /*exit(-1);*/
  /*}*/


  for (i = 0; i < NUM_BUFFERS_PER_THREAD; i++) 
  {
    loops_t.X[i] = ( double* ) malloc ( cf->p * cf->n * cf->x_b * sizeof(double) );
    loops_t.Y[i] = ( double* ) malloc ( cf->NUM_COMPUTE_THREADS * cf->n * sizeof(double) );
    loops_t.B[i] = ( double* ) malloc ( cf->NUM_COMPUTE_THREADS * cf->p * cf->x_b *sizeof(double) );
  }

  loops_t.h = ( double* ) malloc ( cf->t * sizeof(double) );
  Phi = ( double* ) malloc ( cf->n * cf->n * sizeof(double) );
  loops_t.W = ( double* ) malloc ( cf->n * sizeof(double) );
  Z   = ( double* ) malloc ( cf->n * cf->n * sizeof(double) );

  //////// computation specific mem ////////
  loops_t.Winv    = ( double* ) malloc ( cf->NUM_COMPUTE_THREADS * cf->n * sizeof(double) );
  loops_t.XtZWinv = ( double* ) malloc ( cf->NUM_COMPUTE_THREADS * cf->x_b * cf->p * cf->n * sizeof(double)) ;
  loops_t.xtSx    = ( double* ) malloc ( cf->NUM_COMPUTE_THREADS * cf->p * cf->p * sizeof(double) );

  if (loops_t.xtSx == NULL)
  {
	  fprintf(stderr, __FILE__ ": Error, not enough memory\n");
	  exit(EXIT_FAILURE);
  }

  VT_USER_START("PRELOOP");
  Phi_fp = fopen( cf->Phi_path, "rb" );
  read( Phi, Phi_fp, cf->n * cf->n, 0 );
  setenv("OMP_NUM_THREADS", "7", 1);
  preloop(Phi, Z, loops_t.W);
  setenv("OMP_NUM_THREADS", "1", 1);
  free( Phi );
  free( Z );
  fclose( Phi_fp );
  VT_USER_END("PRELOOP");

  loops_t.X_fp = fopen( cf->ZtX_path, "rb");
  loops_t.Y_fp = fopen( cf->ZtY_path, "rb");
  loops_t.B_fp = fopen( cf->B_path, "wb");
  /*if(!x_file) {*/
  /*printf("error opening x_file(%s)! exiting...\n", str_buf);*/
  /*exit(-1);*/
  /*}*/

  /*y_file = fopen(str_buf, "rb");*/
  /*if(!y_file) {*/
  /*printf("error opening y_file(%s)! exiting...\n", str_buf);*/
  /*exit(-1);*/
  /*}*/

  loops_t.sem_io     = (sem_t *) malloc ( cf->NUM_COMPUTE_THREADS * sizeof(sem_t) );
  loops_t.sem_comp_x = (sem_t *) malloc ( cf->NUM_COMPUTE_THREADS * sizeof(sem_t) );
  loops_t.sem_comp_y = (sem_t *) malloc ( cf->NUM_COMPUTE_THREADS * sizeof(sem_t) );
  for (i = 0; i < cf->NUM_COMPUTE_THREADS; i++) 
  {
    sem_init(&loops_t.sem_io[i]    , 0, 0);
    sem_init(&loops_t.sem_comp_x[i], 0, 0);
    sem_init(&loops_t.sem_comp_y[i], 0, 0);
  }

  int value;
  for (i = 0; i < cf->NUM_COMPUTE_THREADS; i++) 
  {
	  sem_getvalue(&loops_t.sem_comp_x[i], &value);
	  /*printf("Initial value of SEM_COMP_X[%d]: %d\n", i, value);*/
	  sem_getvalue(&loops_t.sem_comp_y[i], &value);
	  /*printf("Initial value of SEM_COMP_Y[%d]: %d\n", i, value);*/
      sem_getvalue(&loops_t.sem_io[i], &value);
	  /*printf("Initial value of SEM_IO[i]    : %d\n", value);*/
  }
  
  // only compute by single columns of y
  /*in_p->t_indexed = in_p->t;*/
  /*in_p->y_b = 1;*/
    
  /*printf("starting compute and io threads\n");*/

  VT_USER_START("LOOPS");
  h_fp = fopen( cf->h_path, "r");
  read(loops_t.h, h_fp, cf->t, 0);
  fclose( h_fp );

  loops_t_comp = ( ooc_loops_t* ) malloc ( cf->NUM_COMPUTE_THREADS * sizeof(ooc_loops_t) );
  compute_threads = ( pthread_t * ) malloc ( cf->NUM_COMPUTE_THREADS * sizeof(pthread_t) );
  rc = pthread_create(&io_thread, NULL, io_thread_func, (void*)&loops_t);
  if (rc) {
    printf("error: return code from pthread_create() is %d\n", rc);
    pthread_exit(NULL);
  }

  for (i = 0; i < cf->NUM_COMPUTE_THREADS; i++) {
    memcpy((void*)&loops_t_comp[i], (void*)&loops_t, sizeof(ooc_loops_t));    
    loops_t_comp[i].id = i;

    rc = pthread_create(&compute_threads[i], NULL, compute_thread_func, (void*)&loops_t_comp[i]);
	/*printf("Thread[%d] created\n", i);*/
    if (rc) {
      printf("error: return code from pthread_create() is %d\n", rc);
      pthread_exit(NULL);
    }
  }

  void* retval;
  pthread_join(io_thread, &retval);
  for (i = 0; i < cf->NUM_COMPUTE_THREADS; i++) {
    pthread_join(compute_threads[i], &retval);
  }
  VT_USER_START("END");

  fclose( loops_t.X_fp );
  fclose( loops_t.Y_fp );
  fclose( loops_t.B_fp );

  for (i = 0; i < NUM_BUFFERS_PER_THREAD; i++) {
	  free( loops_t.X[i] );
	  free( loops_t.Y[i] );
	  free( loops_t.B[i] );
  }
  free( loops_t.h );
  free( loops_t.W );

  free( loops_t.Winv );
  free( loops_t.XtZWinv );
  free( loops_t.xtSx );

  for (i = 0; i < cf->NUM_COMPUTE_THREADS; i++) {
	  sem_destroy( &loops_t.sem_comp_x[i] );
	  sem_destroy( &loops_t.sem_comp_y[i] );
      sem_destroy(&loops_t.sem_io[i]);
  }
  free( loops_t.sem_io );
  free( loops_t.sem_comp_x );
  free( loops_t.sem_comp_y );

  free( loops_t_comp );
  free( compute_threads );


  return 0;
}
