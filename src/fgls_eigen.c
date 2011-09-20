#include "fgls.h"
#include "io.h"
#include "mod_x_y.h"

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>

/*#define NUM_COMPUTE_THREADS 1*/
#define NUM_BUFFERS_PER_THREAD 2
/*#define NUM_BUFFERS NUM_BUFFERS_PER_THREAD*NUM_COMPUTE_THREADS*/

double *x[NUM_BUFFERS_PER_THREAD];
double *y[NUM_BUFFERS_PER_THREAD];
double *b[NUM_BUFFERS_PER_THREAD];

double *h;
double *phi, *Z, *W;
double *gWinv, *gXtZWinv, *gxtSx;

/*double *x_compute, *y_compute, *b_compute;*/

sem_t sem_io;
sem_t *sem_comp_x;
sem_t *sem_comp_y;

problem_args* in; //[NUM_COMPUTE_THREADS];

void all_sem_post(sem_t *sem, problem_args *a) {
	int i;
	for ( i = 0; i < a->NUM_COMPUTE_THREADS; i++ )
		sem_post( &sem[i] );
}
void all_sem_wait(sem_t *sem, int index, problem_args *a) {
	int i,
		end = MIN(a->NUM_COMPUTE_THREADS, a->t - index*a->y_b);
	for ( i = 0; i < end; i++ )
		sem_wait( sem );
}
void read_all_x(double *x_cur, int index, problem_args *a) {
	int i;
	for ( i = 0; i < a->NUM_COMPUTE_THREADS; i++ )
	    read_x(&x_cur[i * a->n * a->p * a->x_b], index, a);
}
void read_all_y(double *y_cur, int index, problem_args *a) {
	int i,
		end = MIN(a->NUM_COMPUTE_THREADS, a->t - index*a->y_b);
	for ( i = 0; i < end; i++ )
	    read_y(&y_cur[i * a->n], index+i, a);
}
void write_all_b(double *b, int r, int s, problem_args *a) {
	int i,
		end = MIN(a->NUM_COMPUTE_THREADS, a->t - s*a->y_b);
	for ( i = 0; i < end; i++ )
        write_b(&b[i*a->p*a->x_b], r, s+i, a);
}

void* io_thread_func(void* in) {
  DEF_TIMING();

  double *x_cur  = x[0];
  double *x_next = x[1];
  double *y_cur  = y[0];
  double *y_next = y[1];
  double *b_prev = b[0];
  double *b_cur  = b[1];

  problem_args* args = (problem_args*)in;
  int s, r;
  int i = args->m_indexed;

  BEGIN_TIMING();
  read_all_x(x_cur, 0, args);
  read_all_y(y_cur, 0, args);
  END_TIMING(args->time->io_time);

  all_sem_post(sem_comp_x, args);
  all_sem_post(sem_comp_y, args);
  for (s = 0; s < args->t; s+=args->NUM_COMPUTE_THREADS) {
    swap_buffers(&y_cur, &y_next);
    
    BEGIN_TIMING();
    read_all_y(y_cur, (s + args->NUM_COMPUTE_THREADS) % args->t, args);
    END_TIMING(args->time->io_time);
    all_sem_post(sem_comp_y, args);

    for (r = 0; r < i; r++) {
      swap_buffers(&x_cur, &x_next);

      BEGIN_TIMING();
      read_all_x(x_cur, (r+1) % i, args);
      END_TIMING(args->time->io_time);
      all_sem_post(sem_comp_x, args);
	  /*printf("[IO] Post s: %d r:%d\n", s, r);*/

      BEGIN_TIMING();
      all_sem_wait(&sem_io, s, args);
      END_TIMING(args->time->io_mutex_wait_time);

      swap_buffers(&b_prev, &b_cur);

      BEGIN_TIMING();
      write_all_b(b_prev, r, s, args);
      END_TIMING(args->time->io_time);
    }
  }
  pthread_exit(NULL);
}

void* compute_thread_func(void* in) {
  DEF_TIMING();
  problem_args* args = (problem_args*)in;
  int r, s;
  double *x_cur  = &x[0][args->id * args->p * args->n * args->x_b];
  double *x_next = &x[1][args->id * args->p * args->n * args->x_b];
  double *y_cur  = &y[0][args->id * args->n];
  double *y_next = &y[1][args->id * args->n];
  double *b_prev = &b[0][args->id * args->p * args->x_b];
  double *b_cur  = &b[1][args->id * args->p * args->x_b];

  double *Winv    = &gWinv[args->id * args->n];
  double *XtZWinv = &gXtZWinv[args->id * args->x_b * args->p * args->n];
  double *xtSx    = &gxtSx[args->id * args->p * args->p];

  int i = args->m_indexed;

  int x_inc, y_inc;
  int c, d, l, k;
  int info;

  long temp;
  int mp;
  int pxinc;
  int n = args->n;
  int p = args->p;

  double ONE = 1.0;
  double ZERO = 0.0;
  int iONE = 1;
  int iZERO = 0;

  for (s = args->id; s < args->t; s+=args->NUM_COMPUTE_THREADS) {
    BEGIN_TIMING();
    sem_wait(&sem_comp_y[args->id]);
    END_TIMING(args->time->comp_mutex_wait_time);

    BEGIN_TIMING();
    /* 2) W = sqrt(alpha W - beta I)^-1 */
    // Best order? sqrt - inv
    for (k = 0; k < n; k++)
      Winv[k] = sqrt(1.0 / (h[s]*h[s] * W[k] + (1 - h[s]*h[s])));
    
    /* sqrt(Winv) * ZtY */
    for (l = 0; l < n; l++)
      /*dscal_(&t, &Winv[l], &ZtY[l], &n);*/
      y_cur[l] *= Winv[l];
    END_TIMING(args->time->compute_time);

    for (r = 0; r < i; r++) {
      BEGIN_TIMING();
      sem_wait(&sem_comp_x[args->id]);
      END_TIMING(args->time->comp_mutex_wait_time);

      BEGIN_TIMING();
      x_inc = MIN(args->x_b, args->m - args->x_b*r);
	  mp = p * args->x_b;
      pxinc = p*x_inc;

	  printf("X[0]: %f\n", x_cur[0]);
	  printf("Y[0]: %f\n", y_cur[0]);
      /* X' * Z  * sqrt(Winv) */
      for (l = 0; l < n * mp; l++) XtZWinv[l] = 0.0;
      for (l = 0; l < n; l++)
        daxpy_(&pxinc, &Winv[l], &x_cur[l], &n, &XtZWinv[l*mp], &iONE);
	  printf("W[0]: %f\n", XtZWinv[0]);
      
      /* 7) y = XtZWinv * y */
      dgemv_("N", &pxinc, &n, &ONE, XtZWinv, &mp, y_cur, &iONE, &ZERO, b_cur, &iONE);
	  printf("RHS[0]: %f\n", b_cur[0]);
      
      for (c = 0; c < x_inc; c++)
      {
        /* 5) W = XtZWinv * K^T */
        dsyrk_("L", "N", &p, &n, &ONE, &XtZWinv[c*p], &mp, &ZERO, xtSx, &p);

        /* 8) W^-1 * y */
        dposv_("L", &p, &iONE, xtSx, &p, &b_cur[c*p], &p, &info);
        if (info != 0)
        {
          fprintf(stderr, "Error executing dposv (s: %d, r: %d, c: %d, th: %d): %d\n", s, r, c, args->id, info);
          exit(-1);
        }
      }
      END_TIMING(args->time->compute_time);

      swap_buffers(&x_cur, &x_next);
      swap_buffers(&b_prev, &b_cur);

      sem_post(&sem_io);
    }
    swap_buffers(&y_cur, &y_next);
  }
  pthread_exit(NULL);
}

int fgls_eigen(char* dir, problem_args* in_p) {
  int rc, i;
  pthread_t io_thread;
  pthread_t *compute_threads;

  char str_buf[STR_BUFFER_SIZE];

  sprintf(str_buf, "%s/X.in", dir);
  x_file = fopen(str_buf, "rb");
  if(!x_file) {
    printf("error opening x_file(%s)! exiting...\n", str_buf);
    exit(-1);
  }

  sprintf(str_buf, "%s/X.tmp", dir);
  x_tmp_file = fopen(str_buf, "w+b");
  if(!x_tmp_file) {
    printf("error opening x_tmp_file(%s)! exiting...\n", str_buf);
    exit(-1);
  }

  sprintf(str_buf, "%s/Y.in", dir);
  y_file = fopen(str_buf, "rb");
  if(!y_file) {
    printf("error opening y_file(%s)! exiting...\n", str_buf);
    exit(-1);
  }

  sprintf(str_buf, "%s/Y.tmp", dir);
  y_tmp_file = fopen(str_buf, "w+b");
  if(!y_tmp_file) {
    printf("error opening y_tmp_file(%s)! exiting...\n", str_buf);
    exit(-1);
  }

  sprintf(str_buf, "%s/Phi.in", dir);
  phi_file = fopen(str_buf, "rb");
  if(!phi_file) {
    printf("error opening phi_file(%s)! exiting...\n", str_buf);
    exit(-1);
  }

  sprintf(str_buf, "%s/H.in", dir);
  h_file = fopen(str_buf, "rb");
  if(!h_file) {
    printf("error opening h_file(%s)! exiting...\n", str_buf);
    exit(-1);
  }

  sprintf(str_buf, "%s/B.out", dir);
  b_file = fopen(str_buf, "w+b");
  if(!b_file) {
    printf("error opening b_file(%s)! exiting...\n", str_buf);
    exit(-1);
  }

  for (i = 0; i < NUM_BUFFERS_PER_THREAD; i++) {
    x[i] = (double*)malloc(in_p->NUM_COMPUTE_THREADS * in_p->p * in_p->n * in_p->x_b * sizeof(double));
    y[i] = (double*)malloc(in_p->NUM_COMPUTE_THREADS * in_p->n * sizeof(double));
    b[i] = (double*)malloc(in_p->NUM_COMPUTE_THREADS * in_p->p * in_p->x_b *sizeof(double));
  }
  h   = (double*)malloc(in_p->t * sizeof(double));
  phi = (double*)malloc(in_p->n * in_p->n * sizeof(double));
  W   = (double*)malloc(in_p->n * sizeof(double));
  Z   = (double*)malloc(in_p->n * in_p->n * sizeof(double));

  //////// computation specific mem ////////
  gWinv    = (double*)malloc(in_p->NUM_COMPUTE_THREADS * in_p->n * sizeof(double));
  gXtZWinv = (double*)malloc(in_p->NUM_COMPUTE_THREADS * in_p->x_b * in_p->p * in_p->n * sizeof(double));
  gxtSx    = (double*)malloc(in_p->NUM_COMPUTE_THREADS * in_p->p * in_p->p * sizeof(double));

  if (gxtSx == NULL)
  {
	  fprintf(stderr, "Not enough memory\n");
	  exit(EXIT_FAILURE);
  }
  printf("mem allocated\n");

  read_h(h, 0, in_p);
  read_phi(phi, 0, in_p);
  // (env_variable, value, replace?)
  setenv("OMP_NUM_THREADS", "7", 1);
  eigenDec(in_p->n, phi, Z, W);

  printf("eigen decomp done\n");

  mod_x_y(Z, in_p);
  setenv("OMP_NUM_THREADS", "1", 1);

  printf("modify-x-and-y done\n");

  fclose(x_file);
  fclose(x_tmp_file);
  fclose(y_file);
  fclose(y_tmp_file);

  // use modified X/Y as input
  sprintf(str_buf, "%s/X.tmp", dir);
  x_file = fopen(str_buf, "rb");
  if(!x_file) {
    printf("error opening x_file(%s)! exiting...\n", str_buf);
    exit(-1);
  }

  sprintf(str_buf, "%s/Y.tmp", dir);
  y_file = fopen(str_buf, "rb");
  if(!y_file) {
    printf("error opening y_file(%s)! exiting...\n", str_buf);
    exit(-1);
  }

  sem_init(&sem_io, 0, 0);
  sem_comp_x = (sem_t *) malloc (in_p->NUM_COMPUTE_THREADS * sizeof(sem_t));
  sem_comp_y = (sem_t *) malloc (in_p->NUM_COMPUTE_THREADS * sizeof(sem_t));
  for (i = 0; i < in_p->NUM_COMPUTE_THREADS; i++) {
    sem_init(&sem_comp_x[i], 0, 0);
    sem_init(&sem_comp_y[i], 0, 0);
  }
  
  // only compute by single columns of y
  in_p->t_indexed = in_p->t;
  in_p->y_b = 1;
    
  printf("starting compute and io threads\n");

  in = (problem_args *) malloc (in_p->NUM_COMPUTE_THREADS * sizeof(problem_args));

  rc = pthread_create(&io_thread, NULL, io_thread_func, (void*)&in[0]);
  if (rc) {
    printf("error: return code from pthread_create() is %d\n", rc);
    pthread_exit(NULL);
  }

  compute_threads = (pthread_t *) malloc (in_p->NUM_COMPUTE_THREADS * sizeof(pthread_t));
  for (i = 0; i < in_p->NUM_COMPUTE_THREADS; i++) {
    // copy problem_args, per thread (with thread id)
    memcpy((void*)&in[i], (void*)in_p, sizeof(problem_args));    
    in[i].id = i;
	/*printf("x_b: %d\n", in_p->x_b);*/
	/*printf("i x_b: %d\n", in[i].x_b);*/

    rc = pthread_create(&compute_threads[i], NULL, compute_thread_func, (void*)&in[i]);
    if (rc) {
      printf("error: return code from pthread_create() is %d\n", rc);
      pthread_exit(NULL);
    }
  }

  void* retval;
  pthread_join(io_thread, &retval);
  for (i = 0; i < in_p->NUM_COMPUTE_THREADS; i++) {
    pthread_join(compute_threads[i], &retval);
  }
  printf("joining compute and io threads done\n");

  fclose(x_file);
  fclose(y_file);
  fclose(h_file);
  fclose(phi_file);
  fclose(b_file);

  for (i = 0; i < NUM_BUFFERS_PER_THREAD; i++) {
	  free(x[i]);
	  free(y[i]);
	  free(b[i]);
  }
  free(h);
  free(phi);
  free(W);
  free(Z);

  free(gWinv);
  free(gXtZWinv);
  free(gxtSx);

  free(in);
  free(compute_threads);

  free(sem_comp_x);
  free(sem_comp_y);

  return 0;
}
