#include "fgls.h"
#include "io.h"
#include "mod_x_y.h"

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>

#define NUM_COMPUTE_THREADS 5
#define NUM_BUFFERS_PER_THREAD 2
#define NUM_BUFFERS 2//NUM_BUFFERS_PER_THREAD*NUM_COMPUTE_THREADS

double *x[NUM_BUFFERS];
double *y[NUM_BUFFERS];
double *h[NUM_BUFFERS];
double *b[NUM_BUFFERS];

double *phi, *Z, *W;
double *Winv, *XtZWinv, *xtSx;

double *x_compute, *y_compute, *b_compute;

sem_t sem_io;
sem_t sem_comp_x;
sem_t sem_comp_y;

problem_args in;

void* io_thread_func(void* in) {
  DEF_TIMING();

  double *x_cur = x[0];
  double *x_next = x[1];
  double *y_cur = y[0];
  double *y_next = y[1];
  double *h_cur = h[0];
  double *h_next = h[1];
  double *b_prev = b[0];
  double *b_cur = b[1];

  problem_args* args = (problem_args*)in;
  int s, r;
  int i = args->m_indexed;
  int j = args->t_indexed;

  printf("first io\n");
  BEGIN_TIMING();
  printf("\tx:\n");
  read_x(x_cur, 0, args);
  printf("\ty:\n");
  read_y(y_cur, 0, args);
  printf("\th:\n");
  read_h(h_cur, 0, args);
  END_TIMING(args->time->io_time);

  sem_post(&sem_comp_x);
  sem_post(&sem_comp_y);
  for (s = 0; s < j; s++) {
    swap_buffers(&y_cur, &y_next);
    swap_buffers(&h_cur, &h_next);
    
    BEGIN_TIMING();
    printf("second io\n");
    read_y(y_cur, (s+1) % j, args);
    read_h(h_cur, (s+1) % j, args);
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
  double *h_cur = h[0];
  double *h_next = h[1];
  double *b_prev = b[0];
  double *b_cur = b[1];

  int i = args->m_indexed;
  int j = args->t_indexed;

  int x_inc, y_inc;
  int c, d, l, k;
  int info;

  long temp;
  double *h;
  int mp;
  int n = args->n;
  int p = args->p;

  double ONE = 1.0;
  double ZERO = 0.0;
  int iONE = 1;
  int iZERO = 0;

  printf("compute: start\n");


  for (s = 0; s < j; s++) {
    BEGIN_TIMING();
    sem_wait(&sem_comp_y);
    END_TIMING(args->time->comp_mutex_wait_time);      

    printf("compute: start y_compute\n");
    BEGIN_TIMING();
    y_inc = MIN(args->y_b, args->t - args->y_b*s); 
    h = h_cur;
    double h2;
    for (d = 0; d < y_inc; d++) {
      /* 2) W = sqrt(alpha W - beta I)^-1 */
      // Best order? sqrt - inv
      h2 = h[d]*h[d];
      for (k = 0; k < n; k++)
        Winv[k] = sqrt(1.0 / (h2 * W[k] + (1 - h2)));
      
      /* sqrt(Winv) * ZtY */
      for (l = 0; l < n; l++)
        /*dscal_(&t, &Winv[l], &ZtY[l], &n);*/
        y_cur[(d*n + l)] *= Winv[l];
    }
    END_TIMING(args->time->compute_time);
    for (r = 0; r < i; r++) {
      
      BEGIN_TIMING();
      sem_wait(&sem_comp_x);
      END_TIMING(args->time->comp_mutex_wait_time);      

      BEGIN_TIMING();
      x_inc = MIN(args->x_b, args->m - args->x_b*r);
      mp = p*x_inc;
      {
        /* X' * Z  * sqrt(Winv) */
        for (l = 0; l < n * mp; l++) XtZWinv[l] = 0.0;
        for (l = 0; l < n; l++)
          daxpy_(&mp, &Winv[l], &x_cur[l*mp], &iONE, &XtZWinv[l*mp], &iONE);
        
        /* 7) y = XtZWinv * y */
        dgemv_("N", &mp, &n, &ONE, XtZWinv, &mp, &y_cur[d*n], &iONE, &ZERO, &b_cur[d*mp], &iONE);
        
        for (c = 0; c < x_inc; c++)
        {
          /* 5) W = XtZWinv * K^T */
          dsyrk_("L", "N", &p, &n, &ONE, &XtZWinv[c*p], &mp, &ZERO, xtSx, &p);
          
          /* 8) W^-1 * y */
          dposv_("L", &p, &iONE, xtSx, &p, &b_cur[d*mp + c*p], &p, &info);
          if (info != 0)
          {
            fprintf(stderr, "Error executing dposv: %d\n", info);
            exit(-1);
          }
        }
      }
      END_TIMING(args->time->compute_time);

      swap_buffers(&x_cur, &x_next);
      swap_buffers(&b_prev, &b_cur);

      sem_post(&sem_io);
    }
    swap_buffers(&y_cur, &y_next);
    swap_buffers(&h_cur, &h_next);
  }
  pthread_exit(NULL);
}

int fgls_eigen(char* dir, problem_args* in_p) {
  int rc;
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
  
  memcpy((void*)&in, (void*)in_p, sizeof(problem_args));

  x[0] = (double*)malloc(in.p * in.n * in.x_b * sizeof(double));
  x[1] = (double*)malloc(in.p * in.n * in.x_b * sizeof(double));
  y[0] = (double*)malloc(in.n * in.y_b * sizeof(double));
  y[1] = (double*)malloc(in.n * in.y_b * sizeof(double));
  h[0] = (double*)malloc(in.y_b * sizeof(double));
  h[1] = (double*)malloc(in.y_b * sizeof(double));
  b[0] = (double*)malloc(in.p * in.x_b * in.y_b *sizeof(double));
  b[1] = (double*)malloc(in.p * in.x_b * in.y_b *sizeof(double));

  phi  = (double*)malloc(in.n * in.n * sizeof(double));
  W = (double*) malloc(in.n * sizeof(double));
  Z = (double*) malloc(in.n * in.n * sizeof(double));

  //////// computation specific mem ////////
  //  ZtY = (double *) malloc (n * args->y_b * sizeof(double));
  Winv = (double *) malloc (in.n * sizeof(double));

  //  ZtX = (double *) malloc (args->x_b * p * n * sizeof(double));
  XtZWinv = (double *) malloc (in.x_b * in.p * in.n * sizeof(double));
  xtSx = (double *) malloc (in.p * in.p * sizeof(double));

  printf("mem allocated\n");

  read_phi(phi, 0, &in);
  eigenDec(in.n, phi, Z, W);

  printf("eigen decomp done\n");

  // (env_variable, value, replace?)
  //  setenv("OMP_NUM_THREADS", "7", 1);  
  mod_x_y(Z, &in);
  //  setenv("OMP_NUM_THREADS", "1", 1);  

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

  printf("modify-x-and-y done\n");

  sem_init(&sem_io, 0, 0);
  sem_init(&sem_comp_x, 0, 0);
  sem_init(&sem_comp_y, 0, 0);

  printf("starting compute and io threads\n");

  rc = pthread_create(&io_thread, NULL, io_thread_func, (void*)&in);
  if (rc) {
    printf("error: return code from pthread_create() is %d\n", rc);
    pthread_exit(NULL);
  }
  rc = pthread_create(&compute_thread, NULL, compute_thread_func, (void*)&in);
  if (rc) {
    printf("error: return code from pthread_create() is %d\n", rc);
    pthread_exit(NULL);
  }

  void* retval;
  pthread_join(io_thread, &retval);
  pthread_join(compute_thread, &retval);

  printf("joining compute and io threads\n");

  fclose(x_file);
  fclose(y_file);
  fclose(h_file);
  fclose(phi_file);
  fclose(b_file);

  free(x[0]);
  free(x[1]);
  free(y[0]);
  free(y[1]);

  return 0;
}
