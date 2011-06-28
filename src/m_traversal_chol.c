#include "fgls.h"
#include "io.h"
#include "bio.h"
#if TIMING
#include <sys/time.h>
#include <time.h>
#endif // TIMING

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <unistd.h>
#include <string.h>

double* x[2];
double* y[2];
double* b[2];

double* phi;

double *x_compute, *y_compute, *b_compute;

sem_t sem_io;
sem_t sem_comp;

problem_args in;

void* m_io(void* in) {
#if TIMING
  struct timeval start, end;
#endif // TIMING
  double *x_cur = x[0];
  double *x_next = x[1];
  double *y_cur = y[0];
  double *y_next = y[1];
  double *b_prev = b[0];
  double *b_cur = b[1];
  problem_args* args = (problem_args*)in;
  int s, r;
  int i = args->m_indexed;
  int j = args->m_indexed;
#if DEBUG
  printf("read_x:\n\tx[%d](x=%d)\n", return_buffer_index(x, 2, x_cur), 0);  
  printf("read_y:\n\ty[%d](y=%d)\n", return_buffer_index(y, 2, y_cur), 0);  
#endif // DEBUG
#if TIMING
  gettimeofday(&start, NULL);
#endif // TIMING
  read_x(x_cur, 0, args);
  read_y(y_cur, 0, args);
#if TIMING
  gettimeofday(&end, NULL);
  args->time->io_time += get_diff_ms(&start, &end);
#endif // TIMING
  sem_post(&sem_comp);

  for (s = 0; s < j; s++) {
    swap_buffers(&x_cur, &x_next);

#if DEBUG
    printf("read_x:\n\tx[%d](x=%d)\n", return_buffer_index(x, 2, x_cur), s);        
#endif // DEBUG
#if TIMING
    gettimeofday(&start, NULL);
#endif // TIMING

    read_x(x_cur, (s+1) % j, args);

#if TIMING
    gettimeofday(&end, NULL);
    args->time->io_time += get_diff_ms(&start, &end);
#endif // TIMING

    for (r = 0; r < i; r++) {
      swap_buffers(&y_cur, &y_next);

#if DEBUG
      printf("read_y:\n\ty[%d](y=%d)\n", return_buffer_index(y, 2, y_cur), r);        
#endif // DEBUG
#if TIMING
      gettimeofday(&start, NULL);
#endif // TIMING

      read_y(y_cur, (r+1) % i, args);

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

      swap_buffers(&b_prev, &b_cur);

#if DEBUG
      printf("wrie_b:\n\tb[%d](x=%d)(y=%d)\n", 
             return_buffer_index(b, 2, b_prev), s, r);  
#endif // DEBUG
#if TIMING
      gettimeofday(&start, NULL);
#endif // TIMING

      write_b(b_prev, s, r, args);

#if TIMING
      gettimeofday(&end, NULL);
      args->time->io_time += get_diff_ms(&start, &end);
#endif // TIMING
    }
  }
  pthread_exit(NULL);
}

// compute a (dot) b = c 
//    a - args.p(h) x args.n(k) x args.x_b(i)
//    b - args.n(k) x 1 x args.y_b(j)
//    c - args.p(h) x args.x_b(i) x args.y_b(j)
/*void m_computation( double* a, double* b, double* c, problem_args* args) {
  int i, j, h, k;
  double sum;
  for(j = 0; j < args->y_b; j++) {
    for(i = 0; i < args->x_b; i++) {
      for(h = 0; h < args->p; h++) {
	sum = 0;
	for(k = 0; k < args->n; k++) {
	  sum += ITEM(a, h, args->p, k, args->n, i) *
	    ITEM(b, k, args->n, 0, 1, j);
	}
	ITEM(c, h, args->p, i, args->x_b, j) = sum;
      }
    }
  }
  }*/
void* m_compute(void* in) {
#if TIMING
  struct timeval start, end;
#endif // TIMING
  problem_args* args = (problem_args*)in;
  int s,r;
  double *x_cur = x[0];
  double *x_next = x[1];
  double *y_cur = y[0];
  double *y_next = y[1];
  double *b_prev = b[0];
  double *b_cur = b[1];

  int i = args->t_indexed;
  int j = args->m_indexed;

  for (s = 0; s < j; s++) {
    for (r = 0; r < i; r ++) {
#if TIMING
      gettimeofday(&start, NULL);
#endif // TIMING
      sem_wait(&sem_comp);
#if TIMING
      gettimeofday(&end, NULL);
      args->time->comp_mutex_wait_time += get_diff_ms(&start, &end);
#endif // TIMING
#if DEBUG
      printf("compute(x_cur)\n:");
      print_buffer(x_cur, args->p*args->n*args->x_b);
      printf("compute(y_cur):\n");
      print_buffer(y_cur, args->n*args->y_b);
#endif // DEBUG
#if TIMING
      gettimeofday(&start, NULL);
#endif // TIMING
      bio_chol( args->x_b, args->n, args->p, args->y_b,
                b_prev, x_cur, phi, y_cur,
                args->h );
#if TIMING
      gettimeofday(&end, NULL);
      args->time->compute_time += get_diff_ms(&start, &end);
#endif // TIMING
#if DEBUG
      printf("compute:\n\tx[%d](x=%d)\n\ty[%d](y=%d)\n\tb[%d]\n", 
	     return_buffer_index(x, 2, x_cur), s,
	     return_buffer_index(y, 2, y_cur), r,
	     return_buffer_index(b, 2, b_cur));
      printf("compute(b_cur):\n");
      print_buffer(b_cur, args->p*args->y_b*args->x_b);
#endif // DEBUG

      swap_buffers(&y_cur, &y_next);
      swap_buffers(&b_prev, &b_cur);

      sem_post(&sem_io);
    }
    swap_buffers(&x_cur, &x_next);
  }
  pthread_exit(NULL);
}

int m_traversal_chol(char* x_f, char* y_f, char* phi_f, char* b_f, problem_args* in_p) {
  int rc;
  pthread_t io_thread;
  pthread_t compute_thread;

  x_file = fopen(x_f, "rb");
  if(!x_file) {
    printf("error opening x_file(%s)! exiting...\n", x_f);
    exit(-1);
  }
  y_file = fopen(y_f, "rb");
  if(!y_file) {
    printf("error opening y_file(%s)! exiting...\n", y_f);
    exit(-1);
  }
  phi_file = fopen(phi_f, "rb");
  if(!y_file) {
    printf("error opening phi_file(%s)! exiting...\n", phi_f);
    exit(-1);
  }
  b_file = fopen(b_f, "w+b");
  if(!b_file) {
    printf("error opening b_file(%s)! exiting...\n", b_f);
    exit(-1);
  }
  
  memcpy((void*)&in, (void*)in_p, sizeof(problem_args));

#if DEBUG
  //  write_test_matrices(x_file, y_file, &in);
#endif // DEBUG
  
  x[0] = (double*)malloc(in.p * in.n * in.x_b * sizeof(double));
  x[1] = (double*)malloc(in.p * in.n * in.x_b * sizeof(double));
  y[0] = (double*)malloc(in.n * in.y_b * sizeof(double));
  y[1] = (double*)malloc(in.n * in.y_b * sizeof(double));
  b[0] = (double*)malloc(in.p * in.x_b * in.y_b *sizeof(double));
  b[1] = (double*)malloc(in.p * in.x_b * in.y_b *sizeof(double));
  phi = (double*)malloc(in.n * in.n * sizeof(double));

  sem_init(&sem_io, 0, 0);
  sem_init(&sem_comp, 0, 0);

  read_double(phi, phi_file, in.n, in.n, 0);

  rc = pthread_create(&io_thread, NULL, m_io, (void*)&in);
  if (rc){
    printf("error: return code from pthread_create() is %d\n", rc);
    pthread_exit(NULL);
  }
  rc = pthread_create(&compute_thread, NULL, m_compute, (void*)&in);
  if (rc){
    printf("error: return code from pthread_create() is %d\n", rc);
    pthread_exit(NULL);
  }

  void* retval;

  pthread_join(io_thread, &retval);
  pthread_join(compute_thread, &retval);

#if DEBUG
  print_output(b_file, &in);  
#endif // DEBUG

  fclose(x_file);
  fclose(y_file);
  fclose(phi_file);
  fclose(b_file);

  free(x[0]);
  free(x[1]);
  free(y[0]);
  free(y[1]);
  free(b[0]);
  free(b[1]);
  free(phi);
  return 0;
}
