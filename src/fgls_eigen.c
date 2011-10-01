#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <sys/time.h>
#include <time.h>

#include <pthread.h>
#include <semaphore.h>
#include <aio.h>

#include "blas.h"
#include "lapack.h"
#include "io.h"
#include "timing.h"
#include "fgls_eigen.h"

#if VAMPIR
  #include "vt_user.h"
#endif

#define NUM_BUFFERS_PER_THREAD 2

void swap_buffers(double** b1, double** b2) {
  double* temp;
  temp = *b1;
  *b1 = *b2;
  *b2 = temp;
}

void swap_aiocb(const struct aiocb *x[], const struct aiocb *y[])
{
	const struct aiocb *tmp = x[0];
	x[0] = y[0];
	y[0] = tmp;
}

typedef struct {
	FILE *XR_fp;
	FILE *Y_fp;
	FILE *B_fp;
	double *XL[2];
	double *XL_b;
	double *XLtXL;

	double *X[NUM_BUFFERS_PER_THREAD];
	double *Y[NUM_BUFFERS_PER_THREAD];
	double *B[NUM_BUFFERS_PER_THREAD];

	double *h;
	double *sigma;
	double *W;
	double *Winv;
	double *xtSx;

	FGLS_eigen_t *cf;

	int id;
} ooc_loops_t;


void* compute_thread_func(void* in) 
{
  DEF_TIMING();

  ooc_loops_t *loops_t = ( ooc_loops_t* ) in;
  FGLS_eigen_t *cf = &FGLS_eigen_config;

  int m = cf->m,
	  n = cf->n,
	  p = cf->p,
	  t = cf->t,
	  wXL = cf->wXL,
	  wXR = cf->wXR,
	  x_b = cf->x_b,
	  id = loops_t->id;

  double *xl     = &loops_t->XL[1][id * wXL * n];
  double *xl_b   = &loops_t->XL_b[id * wXL];
  double *xltxl  = &loops_t->XLtXL[id * wXL * wXL];
  double *x_cur  = &loops_t->X[0][id * n * wXR * x_b];
  double *x_next = &loops_t->X[1][id * n * wXR * x_b];
  double *y_cur  = &loops_t->Y[0][id * n];
  double *y_next = &loops_t->Y[1][id * n];
  double *b_cur  = &loops_t->B[0][id * p * x_b];
  double *b_prev = &loops_t->B[1][id * p * x_b];
  double *h      =  loops_t->h;
  double *sigma  =  loops_t->sigma;
  double *W      =  loops_t->W;

  double *Winv = &loops_t->Winv[id * n];
  double *xtSx = &loops_t->xtSx[id * p * p];

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

  struct aiocb aiocb_x_cur,  aiocb_x_next,
			   aiocb_y_cur,  aiocb_y_next,
			   aiocb_b_prev, aiocb_b_cur;
  const struct aiocb *aiocb_x_cur_l[1]  = { &aiocb_x_cur }, 
			         *aiocb_x_next_l[1] = { &aiocb_x_next },
			         *aiocb_y_cur_l[1]  = { &aiocb_y_cur },
			         *aiocb_y_next_l[1] = { &aiocb_y_next },
			         *aiocb_b_prev_l[1] = { &aiocb_b_prev },
			         *aiocb_b_cur_l[1]  = { &aiocb_b_cur };

  bzero( (char *)&aiocb_x_cur,  sizeof(struct aiocb) );
  bzero( (char *)&aiocb_x_next, sizeof(struct aiocb) );
  bzero( (char *)&aiocb_y_cur,  sizeof(struct aiocb) );
  bzero( (char *)&aiocb_y_next, sizeof(struct aiocb) );
  bzero( (char *)&aiocb_b_prev, sizeof(struct aiocb) );
  bzero( (char *)&aiocb_b_cur,  sizeof(struct aiocb) );

  aiocb_x_cur.aio_fildes = fileno( loops_t->XR_fp );
  aiocb_x_cur.aio_buf = x_cur;
  aiocb_x_cur.aio_nbytes = MIN( x_b, m ) * wXR * n * sizeof(double);
  aiocb_x_cur.aio_offset = 0;
  /*aiocb_x_cur.aio_flag = AIO_RAW;*/
  aio_read( &aiocb_x_cur );

  aiocb_y_cur.aio_fildes = fileno( loops_t->Y_fp );
  aiocb_y_cur.aio_buf = y_cur;
  aiocb_y_cur.aio_nbytes = id >= t ? 0 : n * sizeof(double);
  aiocb_y_cur.aio_offset = id >= t ? 0 : id * n * sizeof(double);
  /*aiocb_y_cur.aio_flag = AIO_RAW;*/
  aio_read( &aiocb_y_cur );

  aiocb_x_next.aio_fildes = fileno( loops_t->XR_fp );
  /*aiocb_x_next.aio_flag = AIO_RAW;*/
  aiocb_x_next.aio_buf = x_next;
  aiocb_y_next.aio_fildes = fileno( loops_t->Y_fp );
  /*aiocb_y_next.aio_flag = AIO_RAW;*/
  aiocb_y_next.aio_buf = y_next;

  aiocb_b_prev.aio_fildes = fileno( loops_t->B_fp );
  /*aiocb_b_prev.aio_flag = AIO_RAW;*/
  aiocb_b_prev.aio_buf = b_prev;
  aiocb_b_cur.aio_fildes = fileno( loops_t->B_fp );
  /*aiocb_b_cur.aio_flag = AIO_RAW;*/
  aiocb_b_cur.aio_buf = b_cur;

  int iter = 0;
  struct timeval t0, t1;
  for (j = id; j < t; j += cf->NUM_COMPUTE_THREADS) 
  {
    gettimeofday(&t0, NULL);
    //read( y_cur, loops_t->Y_fp, 
	//		j + cf->NUM_COMPUTE_THREADS > t ? 0 : MIN( cf->NUM_COMPUTE_THREADS * n, (t - (j + cf->NUM_COMPUTE_THREADS)) * n ), 
	//		j + cf->NUM_COMPUTE_THREADS > t ? 0 : (j + cf->NUM_COMPUTE_THREADS) * n );
	struct aiocb *aiocb_y_next_p = (struct aiocb *)aiocb_y_next_l[0];
    bzero( (char *)aiocb_y_next_p,  sizeof(struct aiocb) );

	/*printf("Reading Y = %d in %p\n", 1+j+cf->NUM_COMPUTE_THREADS, y_next);*/
	aiocb_y_next_p->aio_buf = y_next;
	aiocb_y_next_p->aio_fildes = fileno( loops_t->Y_fp );
    aiocb_y_next_p->aio_nbytes = j + cf->NUM_COMPUTE_THREADS >= t ? 0 : n * sizeof(double);
    aiocb_y_next_p->aio_offset = j + cf->NUM_COMPUTE_THREADS >= t ? 0 : (j + cf->NUM_COMPUTE_THREADS) * n * sizeof(double);
	/*aiocb_y.aio_flag = AIO_RAW;*/
	/*printf("Nbytes: %d\n",  aiocb_y_next_p->aio_nbytes);*/
    aio_read( aiocb_y_next_p );

	memcpy( xl, loops_t->XL[0], wXL * n * sizeof(double) );
#if VAMPIR
    VT_USER_START("WAIT_Y");
#endif
	if ( aio_suspend( aiocb_y_cur_l, 1, NULL ) != 0 )
		fprintf(stderr, "Suspend error\n");
	/*printf("Return value: %d\n", aio_return((struct aiocb *)aiocb_y_cur_l[0]));*/
#if VAMPIR
    VT_USER_END("WAIT_Y");
#endif
	/*y_cur = (double *)aiocb_y_cur_l[0]->aio_buf;*/
	/*printf("ZtY(1,   %2d)    = %8f\n", j+1, *y_cur);*/
	/*printf("ZtY(500, %2d)    = %8f\n", j+1, y_cur[499]);*/

#if VAMPIR
    VT_USER_START("COMP_LOOP_Y_CODE");
#endif
  	alpha = sigma[j] * h[j]; // * h[j];
	beta  = sigma[j] * (1 - h[j]); // * h[j]);

    BEGIN_TIMING();
    /* 2) W = sqrt(alpha W - beta I)^-1 */
    // Best order? sqrt - inv
    for (k = 0; k < n; k++)
      Winv[k] = sqrt(1.0 / (alpha * W[k] + beta));
	/*printf("D(1,    1 )     = %8f\n", Winv[0]);*/
	/*printf("D(500, 500)     = %8f\n", Winv[499]);*/
    /* sqrt(Winv) * ZtY */
    for (k = 0; k < n; k++)
      y_cur[k] *= Winv[k];

	/*printf("K(1, %2d)        = %8f\n", j+1, y_cur[0]);*/
	/*printf("K(500, %2d)      = %8f\n", j+1, y_cur[499]);*/

      /* sqrt(Winv) * Z' * XL */
	  int k2;
	  for ( k = 0; k < wXL; k++ )
		  for ( k2 = 0; k2 < n; k2++ )
			  xl[ k * n + k2 ] *= Winv[k2];
	  /*printf("W(1, 1)         = %8f\n", x_cur[0]);*/
	  /*printf("W(500, 3000)    = %8f\n", x_cur[x_b*n*p-1]);*/
      
      /* 7) b = sqrt(Winv) * ZtXL * y */
      dgemv_("T", &n, &wXL, &ONE, xl, &n, y_cur, &iONE, &ZERO, xl_b, &iONE);

    END_TIMING(cf->time->compute_time);
#if VAMPIR
    VT_USER_END("COMP_LOOP_Y_CODE");
#endif
    for (ib = 0; ib < m; ib += x_b) 
	{
#if VAMPIR
		VT_USER_START("READ_X");
#endif
	    struct aiocb *aiocb_x_next_p = (struct aiocb *)aiocb_x_next_l[0];
        bzero( (char *)aiocb_x_next_p,  sizeof(struct aiocb) );

        aiocb_x_next_p->aio_buf = x_next;
	    aiocb_x_next_p->aio_fildes = fileno( loops_t->XR_fp );

		aiocb_x_next_p->aio_nbytes = (ib + x_b) >= m ? MIN( x_b, m ) * wXR * n * sizeof(double) : MIN( x_b * wXR * n, (m - (ib + x_b)) * wXR * n ) * sizeof(double);
		aiocb_x_next_p->aio_offset = (ib + x_b) >= m ? 0 : (ib + x_b) * wXR * n * sizeof(double);
		/*aiocb_x.aio_flag = AIO_RAW;*/
		aio_read( aiocb_x_next_p );
#if VAMPIR
		VT_USER_END("READ_X");
#endif

#if VAMPIR
		VT_USER_START("WAIT_X");
#endif
		BEGIN_TIMING();
		aio_suspend( aiocb_x_cur_l, 1, NULL );
	    END_TIMING(cf->time->comp_mutex_wait_time);
#if VAMPIR
		VT_USER_END("WAIT_X");
#endif
#if VAMPIR
      VT_USER_START("COMP_LOOP_X_CODE");
#endif
      BEGIN_TIMING();
      x_inc = MIN(x_b, m - ib);
	  mpb = p * x_b;
	  mpb_real = p * x_inc;

      /* sqrt(Winv) * Z' * XR */
	  int k2;
	  for ( k = 0; k < x_inc*wXR; k++ )
		  for ( k2 = 0; k2 < n; k2++ )
			  x_cur[ k * n + k2 ] *= Winv[k2];
	  /*printf("W(1, 1)         = %8f\n", x_cur[0]);*/
	  /*printf("W(%d, %d)    = %8f\n", n, x_b*wXR, x_cur[x_b*n*wXR-1]);*/
      
      /* 7) b = WinvZtX * y */
	  /*dgemv_("T", &n, &mpb_real, &ONE, x_cur, &n, y_cur, &iONE, &ZERO, b_cur, &iONE);*/
	  /*printf("RHS(1, %2d)      = %8f\n", j+1, b_cur[0]);*/
	  /*printf("RHS(3000, %2d)   = %8f\n", j+1, b_cur[2999]);*/
      
	  /*dsyrk_("L", "T", &p, &n, &ONE, &x_cur[i*p*n], &n, &ZERO, xtSx, &p);*/
		dsyrk_("L", "T", // LOWER, NO_TRANS, 
				&wXL, &n, // n, k
				&ONE, xl, &n, // KL KL' 
				&ZERO, xltxl, &wXL); // TL of xtSx

      for (i = 0; i < x_inc; i++)
      {
		  int k;
		  memcpy(&b_cur[i*p], xl_b, wXL * sizeof(double));
		  /*printf("DGEMV( T, %d, %d, %2f, %p, %d, %p, %d, %2f, %p, %d);\n",*/
		  /*n, wXR, ONE, &x_cur[i * wXR * n], n, y_cur, iONE, ZERO, b_cur[i*p + wXL], iONE);*/
          dgemv_("T", 
				  &n, &wXR, 
				  &ONE, &x_cur[i * wXR * n], &n, y_cur, &iONE, 
				  &ZERO, &b_cur[i*p + wXL], &iONE);
		  /*printf("RHS(%2d, %2d)    = %8f\n", (ib+i)*p+1,   j+1, b_cur[(ib+i)*p]);*/
		  /*printf("RHS(%2d, %2d)    = %8f\n", (ib+i+1)*p, j+1, b_cur[(ib+i+1)*p-1]);*/
		  /*printf("Thread: %d - Iterating X(%d) - Solving\n", id, i);*/
        /* 5) W = XtZWinv * K^T */
		for( k = 0; k < wXL; k++ )
			dcopy_(&wXL, &xltxl[k*wXL], &iONE, &xtSx[k*p], &iONE);
        dsyrk_("L", "T", 
				&wXR, &n, 
				&ONE, &x_cur[i * wXR * n], &n, 
				&ZERO, &xtSx[wXL * p + wXL], &p);
		dgemm_("T", "N", // NO_TRANS, TRANS, 
				&wXR, &wXL, &n, // m, n, k
				&ONE, &x_cur[i * wXR * n], &n, xl, &n, // KR KL'
				&ZERO, &xtSx[wXL], &p); // xtSx BL

		/*printf("%8f %8f %8f %8f\n", xtSx[0], xtSx[4], xtSx[8], xtSx[12]);*/
		/*printf("%8f %8f %8f %8f\n", xtSx[1], xtSx[5], xtSx[9], xtSx[13]);*/
		/*printf("%8f %8f %8f %8f\n", xtSx[2], xtSx[6], xtSx[10], xtSx[14]);*/
		/*printf("%8f %8f %8f %8f\n", xtSx[3], xtSx[7], xtSx[11], xtSx[15]);*/
        /* 8) W^-1 * y */
        dposv_("L", &p, &iONE, xtSx, &p, &b_cur[i*p], &p, &info);
        if (info != 0)
        {
          fprintf(stderr, "Error executing dposv (y: %d, x: %d, th: %d): %d\n", j, ib+i, id, info);
          exit(-1);
        }
		/*if ( i== 0) printf("res[%d, %2d]: %f\n", ib+i, j, *b_cur);*/
      }
	  /*printf("B(1, %d)         = %8f\n", j+1, b_cur[0]);*/
	  /*printf("B(3000, %d)      = %8f\n", j+1, b_cur[2999]);*/
      END_TIMING(cf->time->compute_time);
#if VAMPIR
      VT_USER_END("COMP_LOOP_X_CODE");
#endif
		
#if VAMPIR
      VT_USER_START("WAIT_B");
#endif
	  if ( iter > 0)
	  {
  		aio_suspend( aiocb_b_prev_l, 1, NULL );
	  }
#if VAMPIR
      VT_USER_END("WAIT_B");
#endif

	  struct aiocb *aiocb_b_cur_p = (struct aiocb *)aiocb_b_cur_l[0];
      bzero( (char *)aiocb_b_cur_p,  sizeof(struct aiocb) );

      aiocb_b_cur_p->aio_buf = b_cur;
	  aiocb_b_cur_p->aio_fildes = fileno( loops_t->B_fp );
	  aiocb_b_cur_p->aio_nbytes = MIN( x_b * p, (m - ib) * p) * sizeof(double);
	  aiocb_b_cur_p->aio_offset = (j * m * p + ib * p) * sizeof(double);
	  /*aiocb_y.aio_flag = AIO_RAW;*/
	  aio_write( aiocb_b_cur_p );

	  swap_aiocb( aiocb_x_cur_l, aiocb_x_next_l );
	  swap_buffers( &x_cur, &x_next);
	  swap_aiocb( aiocb_b_cur_l, aiocb_b_prev_l );
	  swap_buffers( &b_cur, &b_prev);
	  iter++;
    }
    gettimeofday(&t1, NULL);
	printf("Iter time: %ld ms\n", get_diff_ms(&t0, &t1));
	swap_aiocb( aiocb_y_cur_l, aiocb_y_next_l );
	swap_buffers( &y_cur, &y_next);
  }
  aio_suspend( aiocb_x_cur_l, 1, NULL );
  aio_suspend( aiocb_y_cur_l, 1, NULL );
  aio_suspend( aiocb_b_prev_l, 1, NULL );
  pthread_exit(NULL);
}

/*
 * Entry point for the solution of the Feasible Generalized Least-Squares problem
 */
int fgls_eigen( FGLS_eigen_t *cf ) 
{
  int rc, i;
  pthread_t *compute_threads;

  double *Phi,
		 *Z,
		 *W;
  FILE *Phi_fp,
	   *h_fp, *sigma_fp,
	   *XL_fp;
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
    loops_t.X[i] = ( double* ) malloc ( cf->NUM_COMPUTE_THREADS * cf->x_b * cf->p * cf->n * sizeof(double) );
    loops_t.Y[i] = ( double* ) malloc ( cf->NUM_COMPUTE_THREADS * cf->n * sizeof(double) );
    loops_t.B[i] = ( double* ) malloc ( cf->NUM_COMPUTE_THREADS * cf->x_b * cf->p * sizeof(double) );
  }
  // The original in XL[0]
  // The copy to overwrite in XL[1]
  loops_t.XL[0] = ( double* ) malloc ( cf->wXL * cf->n * sizeof(double) );
  loops_t.XL[1] = ( double* ) malloc ( cf->NUM_COMPUTE_THREADS * cf->wXL * cf->n * sizeof(double) );
  loops_t.XL_b  = ( double* ) malloc ( cf->NUM_COMPUTE_THREADS * cf->wXL * sizeof(double) );
  loops_t.XLtXL = ( double* ) malloc ( cf->NUM_COMPUTE_THREADS * cf->wXL * cf->wXL * sizeof(double) );

  loops_t.h = ( double* ) malloc ( cf->t * sizeof(double) );
  loops_t.sigma = ( double* ) malloc ( cf->t * sizeof(double) );
  Phi = ( double* ) malloc ( cf->n * cf->n * sizeof(double) );
  loops_t.W = ( double* ) malloc ( cf->n * sizeof(double) );
  Z   = ( double* ) malloc ( cf->n * cf->n * sizeof(double) );

  loops_t.Winv    = ( double* ) malloc ( cf->NUM_COMPUTE_THREADS * cf->n * sizeof(double) );
  loops_t.xtSx    = ( double* ) malloc ( cf->NUM_COMPUTE_THREADS * cf->p * cf->p * sizeof(double) );

  if (loops_t.xtSx == NULL)
  {
	  fprintf(stderr, __FILE__ ": Error, not enough memory\n");
	  exit(EXIT_FAILURE);
  }

#if VAMPIR
  VT_USER_START("PRELOOP");
#endif
  Phi_fp = fopen( cf->Phi_path, "rb" );
  sync_read( Phi, Phi_fp, cf->n * cf->n, 0 );
  setenv("OMP_NUM_THREADS", "7", 1);
  preloop(Phi, Z, loops_t.W);
  setenv("OMP_NUM_THREADS", "1", 1);
  free( Phi );
  free( Z );
  fclose( Phi_fp );
#if VAMPIR
  VT_USER_END("PRELOOP");
#endif
  loops_t.XR_fp = fopen( cf->ZtXR_path, "rb");
  loops_t.Y_fp = fopen( cf->ZtY_path, "rb");
  loops_t.B_fp = fopen( cf->B_path, "wb");

  XL_fp = fopen( cf->ZtXL_path, "rb" );
  sync_read( loops_t.XL[0], XL_fp, cf->wXL * cf->n, 0 );
  fclose( XL_fp );
  /*if(!x_file) {*/
  /*printf("error opening x_file(%s)! exiting...\n", str_buf);*/
  /*exit(-1);*/
  /*}*/

  /*printf("starting compute and io threads\n");*/

#if VAMPIR
  VT_USER_START("LOOPS");
#endif
  h_fp = fopen( cf->h_path, "r");
  sync_read(loops_t.h, h_fp, cf->t, 0);
  fclose( h_fp );
  sigma_fp = fopen( cf->sigma_path, "r");
  sync_read(loops_t.sigma, sigma_fp, cf->t, 0);
  fclose( sigma_fp );

  printf("LOOPS\n");

  loops_t_comp = ( ooc_loops_t* ) malloc ( cf->NUM_COMPUTE_THREADS * sizeof(ooc_loops_t) );
  compute_threads = ( pthread_t * ) malloc ( cf->NUM_COMPUTE_THREADS * sizeof(pthread_t) );
  for (i = 0; i < cf->NUM_COMPUTE_THREADS; i++) 
  {
    memcpy((void*)&loops_t_comp[i], (void*)&loops_t, sizeof(ooc_loops_t));    
    loops_t_comp[i].id = i;

    rc = pthread_create(&compute_threads[i], NULL, compute_thread_func, (void*)&loops_t_comp[i]);
    if (rc) 
	{
      fprintf(stderr, "ERROR in " __FILE__ ": pthread_create() returned %d\n", rc);
      pthread_exit(NULL);
    }
  }

  void* retval;
  for (i = 0; i < cf->NUM_COMPUTE_THREADS; i++)
    pthread_join(compute_threads[i], &retval);
#if VAMPIR
  VT_USER_START("END");
#endif

  printf("LOOPS Done\n");

  fclose( loops_t.XR_fp );
  fclose( loops_t.Y_fp );
  fclose( loops_t.B_fp );

  for (i = 0; i < NUM_BUFFERS_PER_THREAD; i++) 
  {
	  free( loops_t.X[i] );
	  free( loops_t.Y[i] );
	  free( loops_t.B[i] );
  }
  free( loops_t.XL[0] );
  free( loops_t.XL[1] );
  free( loops_t.XL_b  );
  free( loops_t.XLtXL );

  free( loops_t.h );
  free( loops_t.sigma );
  free( loops_t.W );

  free( loops_t.Winv );
  free( loops_t.xtSx );

  free( loops_t_comp );
  free( compute_threads );


  return 0;
}
