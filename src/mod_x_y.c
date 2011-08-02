#include "fgls.h"
#include "io.h"

#include "options.h"

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>

double* x[2];
double* y[2];

double* Z;

double *x_compute, *y_compute, *b_compute;

sem_t sem_io;
sem_t sem_comp;

problem_args in;

void* io_x(void* in) {
  DEF_TIMING();

  double *x_cur = x[0];
  double *x_next = x[1];
  problem_args* args = (problem_args*)in;
  int r;
  int i = args->m_indexed;
  int j = args->t_indexed;

  BEGIN_TIMING();
  read_x(x_cur, 0, args);
  END_TIMING(args->time->io_time);

  sem_post(&sem_comp);

  for (r = 0; r < i; r++) {
    swap_buffers(&x_cur, &x_next);
    
    BEGIN_TIMING();
    read_x(x_cur, (r+1) % i, args);
    END_TIMING(args->time->io_time);

    sem_post(&sem_comp);
    
    BEGIN_TIMING();
    sem_wait(&sem_io);
    END_TIMING(args->time->io_mutex_wait_time);

    BEGIN_TIMING();
    write_x(x_cur, r, args);
    END_TIMING(args->time->io_time);

  }
  pthread_exit(NULL);
}
void* io_y(void* in) {
  DEF_TIMING();

  double *y_cur = y[0];
  double *y_next = y[1];
  problem_args* args = (problem_args*)in;
  int s;
  int j = args->t_indexed;

  BEGIN_TIMING();
  read_y(y_cur, 0, args);
  END_TIMING(args->time->io_time);

  sem_post(&sem_comp);

  for (s = 0; s < j; s++) {
    swap_buffers(&y_cur, &y_next);
    
    BEGIN_TIMING();
    read_y(y_cur, (s+1) % j, args);
    END_TIMING(args->time->io_time);
    
    sem_post(&sem_comp);
    
    BEGIN_TIMING();
    sem_wait(&sem_io);
    END_TIMING(args->time->io_mutex_wait_time);
    
    BEGIN_TIMING();
    write_y(y_cur, s, args);
    END_TIMING(args->time->io_time);
  }
  pthread_exit(NULL);
}

void* compute_x(void* in) {
  DEF_TIMING();

  problem_args* args = (problem_args*)in;
  int r;
  double *x_cur = x[0];
  double *x_next = x[1];

  double ONE = 1.0;
  double ZERO = 0.0;

  int i = args->m_indexed;
  int j = args->t_indexed;

  int mp, xp, y_b;
  long temp;
  for (r = 0; r < i; r++) {
    BEGIN_TIMING();
    sem_wait(&sem_comp);
    END_TIMING(args->time->comp_mutex_wait_time);

    BEGIN_TIMING();
    mp = (args->x_b*args->p);
    xp = MIN(args->x_b, args->m - args->x_b*r)*args->p;
    /* X <- X^T * Z */
    dgemm_("T", "N", &xp, &args->n, &args->n, &ONE, x_cur, &args->n, Z, &args->n, &ZERO, x_cur, &mp);
    END_TIMING(args->time->compute_time);

    swap_buffers(&x_cur, &x_next);
  }
  pthread_exit(NULL);
}

void* compute_y(void* in) {
  DEF_TIMING();

  problem_args* args = (problem_args*)in;
  int s;
  double *y_cur = y[0];
  double *y_next = y[1];

  double ONE = 1.0;
  double ZERO = 0.0;

  int j = args->t_indexed;

  int mp, xp, y_b;
  long temp;
  for (s = 0; s < j; s++) {
    BEGIN_TIMING();
    sem_wait(&sem_comp);
    END_TIMING(args->time->comp_mutex_wait_time);

    BEGIN_TIMING();
    y_b = MIN(args->y_b, args->t - args->y_b * s);
    /* 6) ZtY = Z' * Y */
    dgemm_("T", "N", &args->n, &y_b, &args->n, &ONE, Z, &args->n, y_cur, &args->n, &ZERO, y_cur, &args->n);
    END_TIMING(args->time->compute_time);

    swap_buffers(&y_cur, &y_next);
  }
  pthread_exit(NULL);
}

int mod_x_y(double* Z_in, problem_args* in_p) {
  int rc;
  pthread_t io_thread;
  pthread_t compute_thread;
  problem_args in;

  memcpy((void*)&in, (void*)in_p, sizeof(problem_args));

  x[0] = (double*)malloc(in.p * in.n * in.x_b * sizeof(double));
  x[1] = (double*)malloc(in.p * in.n * in.x_b * sizeof(double));
  y[0] = (double*)malloc(in.n * in.y_b * sizeof(double));
  y[1] = (double*)malloc(in.n * in.y_b * sizeof(double));

  Z = Z_in;

  sem_init(&sem_io, 0, 0);
  sem_init(&sem_comp, 0, 0);

  rc = pthread_create(&io_thread, NULL, io_x, (void*)&in);
  if (rc){
    printf("error: return code from pthread_create() is %d\n", rc);
    pthread_exit(NULL);
  }
  rc = pthread_create(&compute_thread, NULL, compute_x, (void*)&in);
  if (rc){
    printf("error: return code from pthread_create() is %d\n", rc);
    pthread_exit(NULL);
  }

  void* retval;

  pthread_join(io_thread, &retval);
  pthread_join(compute_thread, &retval);

  // reset semaphores
  sem_init(&sem_io, 0, 0);
  sem_init(&sem_comp, 0, 0);

  rc = pthread_create(&io_thread, NULL, io_y, (void*)&in);
  if (rc){
    printf("error: return code from pthread_create() is %d\n", rc);
    pthread_exit(NULL);
  }
  rc = pthread_create(&compute_thread, NULL, compute_y, (void*)&in);
  if (rc){
    printf("error: return code from pthread_create() is %d\n", rc);
    pthread_exit(NULL);
  }

  pthread_join(io_thread, &retval);
  pthread_join(compute_thread, &retval);

  free(x[0]);
  free(x[1]);
  free(y[0]);
  free(y[1]);

  return 0;
}
