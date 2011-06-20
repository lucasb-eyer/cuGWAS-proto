#include "fgls.h"
#include "io.h"

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <unistd.h>

double* x[2];
double* y[2];
double* b[2];

double *x_compute, *y_compute, *b_compute;

#if USE_NAMED_SEMAPHORES
sem_t* sem_io;
sem_t* sem_comp;
#else
sem_t sem_io;
sem_t sem_comp;
#endif

#if USE_NAMED_SEMAPHORES
#define THREAD_POOL_SEMAPHORE_NAME "fgls-thread-pool"
#define FORKER_SEMAPHORE_NAME_TEMPLATE "fgls_thread_%d"
#define FORKER_SEMAPHORE_NAME_BYTES_NEEDED 16/* Bytes */
#define NAMED_SEMAPHORE_PERMISSIONS 0777
#endif

problem_args in;
/*
#ifdef DEBUG
  printf("read_y:\n\ty[%d](y=%d)\n", return_buffer_index(y, 2, buf), index);  
  #endif*/
void* io(void* in) {
  double *x_cur = x[0];
  double *x_next = x[1];
  double *y_cur = y[0];
  double *y_next = y[1];
  double *b_prev = b[0];
  double *b_cur = b[1];
  problem_args* args = (problem_args*)in;
  int s, r;
  int i = args->t_indexed;
  int j = args->m_indexed;
#ifdef DEBUG
  printf("read_x:\n\tx[%d](x=%d)\n", return_buffer_index(x, 2, x_cur), 0);  
  printf("read_y:\n\ty[%d](y=%d)\n", return_buffer_index(y, 2, y_cur), 0);  
#endif
  read_x(x_cur, 0, args);
  read_y(y_cur, 0, args);

#if USE_NAMED_SEMAPHORES
  sem_post(sem_comp);
#else
  sem_post(&sem_comp);
#endif

  for (r = 0; r < i; r++) {
    swap_buffers(&y_cur, &y_next);
#ifdef DEBUG
  printf("read_y:\n\ty[%d](y=%d)\n", return_buffer_index(y, 2, y_cur), r);  
#endif
    read_y(y_cur, (r+1) % i, args);
    for (s = 0; s < j; s++) {
      swap_buffers(&x_cur, &x_next);
#ifdef DEBUG
  printf("read_x:\n\tx[%d](x=%d)\n", return_buffer_index(x, 2, x_cur), s);  
#endif
      read_x(x_cur, (s+1) % j, args);
#if USE_NAMED_SEMAPHORES
      sem_post(sem_comp);
      sem_wait(sem_io);
#else
      sem_post(&sem_comp);
      sem_wait(&sem_io);
#endif
      swap_buffers(&b_prev, &b_cur);
#ifdef DEBUG
      printf("write_b:\n\tb[%d](x=%d)(y=%d)\n", 
             return_buffer_index(b, 2, b_prev), s, r);  
#endif
      write_b(b_prev, s, r, args);
    }
  }
  pthread_exit(NULL);
}

// compute a (dot) b = c 
//    a - args.p(h) x args.n(k) x args.x_b(i)
//    b - args.n(k) x 1 x args.y_b(j)
//    c - args.p(h) x args.x_b(i) x args.y_b(j)
void computation( double* a, double* b, double* c, problem_args* args) {
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
}
void* compute(void* in) {
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

  for (r = 0; r < i; r++) {
    for (s = 0; s < j; s ++) {

#if USE_NAMED_SEMAPHORES
      sem_wait(sem_comp);
#else
      sem_wait(&sem_comp);
#endif

#ifdef DEBUG
      printf("compute(x_cur)\n:");
      print_buffer(x_cur, args->p*args->n*args->x_b);
      printf("compute(y_cur):\n");
      print_buffer(y_cur, args->n*args->y_b);
#endif

      computation(x_cur, y_cur, b_cur, args);
#ifdef DEBUG
      printf("compute:\n\tx[%d](x=%d)\n\ty[%d](y=%d)\n\tb[%d]\n", 
	     return_buffer_index(x, 2, x_cur), s,
	     return_buffer_index(y, 2, y_cur), r,
	     return_buffer_index(b, 2, b_cur));
      printf("compute(b_cur):\n");
      print_buffer(b_cur, args->p*args->y_b*args->x_b);
#endif

      swap_buffers(&x_cur, &x_next);
      swap_buffers(&b_prev, &b_cur);

#if USE_NAMED_SEMAPHORES
      sem_post(sem_io);
#else
      sem_post(&sem_io);
#endif

    }
    swap_buffers(&y_cur, &y_next);
  }
  pthread_exit(NULL);

}

int t_traversal(int argc, char* argv[]) {
  int rc;
  pthread_t io_thread;
  pthread_t compute_thread;

  if (argc != 4) {
    printf("usage: %s <x-in-file> <y-in-file> <b-out-file>\n", argv[0]);
    exit(-1);
  }
  x_file = fopen(argv[1], "w+b");
  if(!x_file) {
    printf("error opening x_file(%s)! exiting...\n", argv[1]);
    exit(-1);
  }
  y_file = fopen(argv[2], "w+b");
  if(!y_file) {
    printf("error opening y_file(%s)! exiting...\n", argv[2]);
    exit(-1);
  }
  b_file = fopen(argv[1], "w+b");
  if(!b_file) {
    printf("error opening b_file(%s)! exiting...\n", argv[3]);
    exit(-1);
  }

  printf("Please enter parameters\n");
  printf("\tm: ");
  scanf("%d", &in.m);
  printf("\tt: ");
  scanf("%d", &in.t);
  printf("\tm blocksize: ");
  scanf("%d", &in.x_b);
  printf("\tt blocksize: ");
  scanf("%d", &in.y_b);
  printf("\tn: ");
  scanf("%d", &in.n);
  printf("\tp: ");
  scanf("%d", &in.p);

  in.m_indexed = (int) ((double) in.m/in.x_b);
  in.t_indexed = (int) ((double) in.t/in.y_b);
  
  write_test_matrices(x_file, y_file, &in);
  
  x[0] =  (double*)malloc(in.p * in.n * in.x_b * sizeof(double));
  x[1] = (double*)malloc(in.p * in.n * in.x_b * sizeof(double));
  y[0] =  (double*)malloc(in.n * in.y_b * sizeof(double));
  y[1] = (double*)malloc(in.n * in.y_b * sizeof(double));
  b[0] = (double*)malloc(in.p * in.x_b * in.y_b *sizeof(double));
  b[1] =  (double*)malloc(in.p * in.x_b * in.y_b *sizeof(double));
  
#if USE_NAMED_SEMAPHORES
  sem_io = sem_open(THREAD_POOL_SEMAPHORE_NAME, O_CREAT, 
                    NAMED_SEMAPHORE_PERMISSIONS, 0);
  sem_comp = sem_open(THREAD_POOL_SEMAPHORE_NAME, O_CREAT, 
                      NAMED_SEMAPHORE_PERMISSIONS, 0);
#else
  sem_init(&sem_io, 0, 0);
  sem_init(&sem_comp, 0, 0);
#endif

  rc = pthread_create(&io_thread, NULL, io, (void*)&in);
  if (rc){
    printf("error: return code from pthread_create() is %d\n", rc);
    pthread_exit(NULL);
  }
  rc = pthread_create(&compute_thread, NULL, compute, (void*)&in);
  if (rc){
    printf("error: return code from pthread_create() is %d\n", rc);
    pthread_exit(NULL);
  }

  void* retval;

  pthread_join(io_thread, &retval);
  pthread_join(compute_thread, &retval);

  fclose(x_file);
  fclose(y_file);
  fclose(b_file);

  print_output(&in);

  free(x[0]);
  free(x[1]);
  free(y[0]);
  free(y[1]);
  free(b[0]);
  free(b[1]);

  pthread_exit(NULL);
  return 0;
}
