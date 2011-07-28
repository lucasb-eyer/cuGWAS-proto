#include "fgls.h"
#include "io.h"

#if TIMING
#include <sys/time.h>
#include <time.h>
#endif // TIMING

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>

#include <string.h>

double* x[2];
double* y[2];
double* h[2];
double* b[2];

double* phi, *Z, *W;

double *x_compute, *y_compute, *b_compute;

sem_t sem_io;
sem_t sem_comp;

problem_args in;

void* io_x_y(void* in) {
#if TIMING
  struct timeval start, end;
#endif // TIMING
  double *x_cur = x[0];
  double *x_next = x[1];
  double *y_cur = y[0];
  double *y_next = y[1];
  problem_args* args = (problem_args*)in;
  int s, r;
  int i = args->m_indexed;
  int j = args->t_indexed;

#if TIMING
  gettimeofday(&start, NULL);
#endif // TIMING
  read_x(x_cur, 0, args);
#if TIMING
  gettimeofday(&end, NULL);
  args->time->io_time += get_diff_ms(&start, &end);
#endif // TIMING
  sem_post(&sem_comp);

  for (r = 0; r < i; r++) {
    swap_buffers(&x_cur, &x_next);
    
#if TIMING
    gettimeofday(&start, NULL);
#endif // TIMING
    
    read_x(x_cur, (r+1) % i, args);
    
#if TIMING
    gettimeofday(&end, NULL);
    args->time->io_time += get_diff_ms(&start, &end);
#endif // TIMING
    
    sem_post(&sem_comp);
    
#if TIMING
    gettimeofday(&start, NULL);
#endif // TIMING
    
    sem_wait(&sem_io);
    
#if TIMING
    gettimeofday(&end, NULL);
    args->time->io_mutex_wait_time += get_diff_ms(&start, &end);
#endif // TIMING
    
#if TIMING
    gettimeofday(&start, NULL);
#endif // TIMING
    
    write_x(x_cur, r, args);
    
#if TIMING
    gettimeofday(&end, NULL);
    args->time->io_time += get_diff_ms(&start, &end);
#endif // TIMING
    
  }

#if TIMING
  gettimeofday(&start, NULL);
#endif // TIMING
  read_y(y_cur, 0, args);
#if TIMING
  gettimeofday(&end, NULL);
  args->time->io_time += get_diff_ms(&start, &end);
#endif // TIMING
  sem_post(&sem_comp);

  for (s = 0; s < j; s++) {
    swap_buffers(&y_cur, &y_next);
    
#if TIMING
    gettimeofday(&start, NULL);
#endif // TIMING
    
    read_y(y_cur, (s+1) % j, args);
    
#if TIMING
    gettimeofday(&end, NULL);
    args->time->io_time += get_diff_ms(&start, &end);
#endif // TIMING
    
    sem_post(&sem_comp);
    
#if TIMING
    gettimeofday(&start, NULL);
#endif // TIMING
    
    sem_wait(&sem_io);
    
#if TIMING
    gettimeofday(&end, NULL);
    args->time->io_mutex_wait_time += get_diff_ms(&start, &end);
#endif // TIMING
    
#if TIMING
    gettimeofday(&start, NULL);
#endif // TIMING
    
    write_y(y_cur, r, args);
    
#if TIMING
    gettimeofday(&end, NULL);
    args->time->io_time += get_diff_ms(&start, &end);
#endif // TIMING
    
  }
  pthread_exit(NULL);
}

void* compute_x_y(void* in) {
#if TIMING
  struct timeval start, end;
#endif // TIMING
  problem_args* args = (problem_args*)in;
  int r, s;
  double *x_cur = x[0];
  double *x_next = x[1];
  double *y_cur = y[0];
  double *y_next = y[1];

  int i = args->m_indexed;
  int j = args->t_indexed;

  int x_inc, y_inc;
  long temp;
  for (r = 0; r < i; r++) {
#if TIMING
      gettimeofday(&start, NULL);
#endif // TIMING
      sem_wait(&sem_comp);
#if TIMING
      gettimeofday(&end, NULL);
      temp =  get_diff_ms(&start, &end);
      args->time->comp_mutex_wait_time += temp;
#endif // TIMING

#if TIMING
      gettimeofday(&start, NULL);
#endif // TIMING
      //      x_inc = MIN(args->x_b, args->m - args->x_b*r);
      //      y_inc = MIN(args->y_b, args->t - args->y_b*s); 
      //      bio_eigen(x_inc, args->n, args->p, y_inc,
      //                b_cur, x_cur, Z, W, y_cur,
      //                h_cur);
#if TIMING
      gettimeofday(&end, NULL);
      args->time->compute_time += get_diff_ms(&start, &end);
#endif // TIMING

    swap_buffers(&x_cur, &x_next);
  }

  for (s = 0; s < j; s++) {
#if TIMING
      gettimeofday(&start, NULL);
#endif // TIMING
      sem_wait(&sem_comp);
#if TIMING
      gettimeofday(&end, NULL);
      temp =  get_diff_ms(&start, &end);
      args->time->comp_mutex_wait_time += temp;
#endif // TIMING

#if TIMING
      gettimeofday(&start, NULL);
#endif // TIMING
      //      x_inc = MIN(args->x_b, args->m - args->x_b*r);
      //      y_inc = MIN(args->y_b, args->t - args->y_b*s); 
      //      bio_eigen(x_inc, args->n, args->p, y_inc,
      //                b_cur, x_cur, Z, W, y_cur,
      //                h_cur);
#if TIMING
      gettimeofday(&end, NULL);
      args->time->compute_time += get_diff_ms(&start, &end);
#endif // TIMING

    swap_buffers(&y_cur, &y_next);
  }
  pthread_exit(NULL);
}

int mod_x_y(double* Z, problem_args* in_p) {
  int rc;
  pthread_t io_thread;
  pthread_t compute_thread;
  problem_args in;

  memcpy((void*)&in, (void*)in_p, sizeof(problem_args));

  x[0] = (double*)malloc(in.p * in.n * in.x_b * sizeof(double));
  x[1] = (double*)malloc(in.p * in.n * in.x_b * sizeof(double));
  y[0] = (double*)malloc(in.n * in.y_b * sizeof(double));
  y[1] = (double*)malloc(in.n * in.y_b * sizeof(double));

  sem_init(&sem_io, 0, 0);
  sem_init(&sem_comp, 0, 0);

  rc = pthread_create(&io_thread, NULL, io_x_y, (void*)&in);
  if (rc){
    printf("error: return code from pthread_create() is %d\n", rc);
    pthread_exit(NULL);
  }
  rc = pthread_create(&compute_thread, NULL, compute_x_y, (void*)&in);
  if (rc){
    printf("error: return code from pthread_create() is %d\n", rc);
    pthread_exit(NULL);
  }

  void* retval;

  pthread_join(io_thread, &retval);
  pthread_join(compute_thread, &retval);

  free(x[0]);
  free(x[1]);
  free(y[0]);
  free(y[1]);

  return 0;
}
