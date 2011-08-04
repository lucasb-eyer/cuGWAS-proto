#include "fgls.h"
#include "io.h"
#include "mod_x_y.h"

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>

#define NUM_COMPUTE_THREADS 1
#define NUM_BUFFERS_PER_THREAD 2
#define NUM_BUFFERS NUM_BUFFERS_PER_THREAD*NUM_COMPUTE_THREADS

double *x[NUM_BUFFERS];
double *y[NUM_BUFFERS];
double *b[NUM_BUFFERS];

double *h;
double *phi, *Z, *W;
double *Winv, *XtZWinv, *xtSx;

double *x_compute, *y_compute, *b_compute;

sem_t sem_io;
sem_t sem_comp_x;
sem_t sem_comp_y;

problem_args in[NUM_COMPUTE_THREADS];

void* io_thread_func(void* in) {
  DEF_TIMING();

  double *x_cur = x[0];
  double *x_next = x[1];
  double *y_cur = y[0];
  double *y_next = y[1];
  double *b_prev = b[0];
  double *b_cur = b[1];

  problem_args* args = (problem_args*)in;
  int s, r;
  int i = args->m_indexed;

  printf("first io\n");
  BEGIN_TIMING();
  printf("x\n");
  read_x(x_cur, 0, args);
  printf("y\n");
  read_y(y_cur, 0, args);
  END_TIMING(args->time->io_time);

  sem_post(&sem_comp_x);
  sem_post(&sem_comp_y);
  for (s = 0; s < args->t; s++) {
    swap_buffers(&y_cur, &y_next);
    
    BEGIN_TIMING();
    read_y(y_cur, (s+1) % args->t, args);
    END_TIMING(args->time->io_time);
    sem_post(&sem_comp_y);

    for (r = 0; r < i; r++) {
      swap_buffers(&x_cur, &x_next);

      BEGIN_TIMING();
      read_x(x_cur, (r+1) % i, args);
      END_TIMING(args->time->io_time);
      sem_post(&sem_comp_x);

      BEGIN_TIMING();
      sem_wait(&sem_io);
      END_TIMING(args->time->io_mutex_wait_time);

      swap_buffers(&b_prev, &b_cur);

      BEGIN_TIMING();
      write_b(b_prev, r, s, args);
      END_TIMING(args->time->io_time);
    }
  }
  pthread_exit(NULL);
}

void* compute_thread_func(void* in) {
  DEF_TIMING();
  problem_args* args = (problem_args*)in;
  int r, s;
  double *x_cur = x[0];
  double *x_next = x[1];
  double *y_cur = y[0];
  double *y_next = y[1];
  double *b_prev = b[0];
  double *b_cur = b[1];

  int i = args->m_indexed;

  int x_inc, y_inc;
  int c, d, l, k;
  int info;

  long temp;
  int mp;
  int n = args->n;
  int p = args->p;

  double ONE = 1.0;
  double ZERO = 0.0;
  int iONE = 1;
  int iZERO = 0;

  printf("compute: start\n");

  for (s = 0; s < args->t; s++) {
    BEGIN_TIMING();
    sem_wait(&sem_comp_y);
    END_TIMING(args->time->comp_mutex_wait_time);

    BEGIN_TIMING();
    /* 2) W = sqrt(alpha W - beta I)^-1 */
    // Best order? sqrt - inv
    for (k = 0; k < n; k++)
      Winv[k] = sqrt(1.0 / (h[s]*h[s] * W[k] + (1 - h[s]*h[s])));
    
    /* sqrt(Winv) * ZtY */
    for (l = 0; l < n; l++)
      /*dscal_(&t, &Winv[l], &ZtY[l], &n);*/
      y_cur[s*n + l] *= Winv[l];
    END_TIMING(args->time->compute_time);

    for (r = 0; r < i; r++) {
      BEGIN_TIMING();
      sem_wait(&sem_comp_x);
      END_TIMING(args->time->comp_mutex_wait_time);

      BEGIN_TIMING();
      x_inc = MIN(args->x_b, args->m - args->x_b*r);
      mp = p*x_inc;

      /* X' * Z  * sqrt(Winv) */
      for (l = 0; l < n * mp; l++) XtZWinv[l] = 0.0;
      for (l = 0; l < n; l++)
        daxpy_(&mp, &Winv[l], &x_cur[l*mp], &iONE, &XtZWinv[l*mp], &iONE);
      
      /* 7) y = XtZWinv * y */
      dgemv_("N", &mp, &n, &ONE, XtZWinv, &mp, &y_cur[s*n], &iONE, &ZERO, &b_cur[s*mp], &iONE);
      
      for (c = 0; c < x_inc; c++)
      {
        /* 5) W = XtZWinv * K^T */
        dsyrk_("L", "N", &p, &n, &ONE, &XtZWinv[c*p], &mp, &ZERO, xtSx, &p);

        printf("x_cur:\n");
        print_buffer(x_cur, mp);
        printf("y_cur:\n");
        print_buffer(y_cur, n);
        printf("xtSx:\n");
        print_buffer(xtSx, p);
        printf("b_cur:\n");
        print_buffer(&b_cur[s*mp + c*p], p);
        /* 8) W^-1 * y */
        dposv_("L", &p, &iONE, xtSx, &p, &b_cur[s*mp + c*p], &p, &info);
        if (info != 0)
        {
          fprintf(stderr, "Error executing dposv: %d\n", info);
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
  pthread_t compute_thread;

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

  for (i = 0; i < NUM_BUFFERS; i++) {
    x[i] = (double*)malloc(in_p->p * in_p->n * in_p->x_b * sizeof(double));
    y[i] = (double*)malloc(in_p->n * sizeof(double));
    b[i] = (double*)malloc(in_p->p * in_p->x_b *sizeof(double));
  }
  h   = (double*)malloc(in_p->t * sizeof(double));
  phi = (double*)malloc(in_p->n * in_p->n * sizeof(double));
  W   = (double*)malloc(in_p->n * sizeof(double));
  Z   = (double*)malloc(in_p->n * in_p->n * sizeof(double));

  //////// computation specific mem ////////
  Winv    = (double*)malloc(in_p->n * sizeof(double));
  XtZWinv = (double*)malloc(in_p->x_b * in_p->p * in_p->n * sizeof(double));
  xtSx    = (double*)malloc(in_p->p * in_p->p * sizeof(double));

  printf("mem allocated\n");

  read_h(h, 0, in_p);
  read_phi(phi, 0, in_p);
  eigenDec(in_p->n, phi, Z, W);

  printf("eigen decomp done\n");

  // (env_variable, value, replace?)
  setenv("OMP_NUM_THREADS", "7", 1);
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

  // for (i = 0; i < NUM_COMPUTE_THREADS; i++) {
  sem_init(&sem_io, 0, 0);
  sem_init(&sem_comp_x, 0, 0);
  sem_init(&sem_comp_y, 0, 0);
    //  }
  
  // only compute by single columns of y
  in_p->t_indexed = in_p->t;
  in_p->y_b = 1;
    
  printf("starting compute and io threads\n");

  rc = pthread_create(&io_thread, NULL, io_thread_func, (void*)&in);
  if (rc) {
    printf("error: return code from pthread_create() is %d\n", rc);
    pthread_exit(NULL);
  }

  for (i = 0; i < NUM_COMPUTE_THREADS; i++) {
    // copy problem_args, per thread (with thread id)
    memcpy((void*)&in[i], (void*)in_p, sizeof(problem_args));    
    in[i].id = i;

    rc = pthread_create(&compute_thread, NULL, compute_thread_func, (void*)&in);
    if (rc) {
      printf("error: return code from pthread_create() is %d\n", rc);
      pthread_exit(NULL);
    }
  }

  void* retval;
  pthread_join(io_thread, &retval);
  for (i = 0; i < NUM_COMPUTE_THREADS; i++) {
    pthread_join(compute_thread, &retval);
  }
  printf("joining compute and io threads done\n");

  fclose(x_file);
  fclose(y_file);
  fclose(h_file);
  fclose(phi_file);
  fclose(b_file);

  for (i = 0; i < NUM_BUFFERS; i ++) {
    free(x[i]);
    free(y[i]);
    free(b[i]);
  }
  free(h);
  free(phi);
  free(W);
  free(Z);

  free(Winv);
  free(XtZWinv);
  free(xtSx);

  return 0;
}
