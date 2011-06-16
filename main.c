#define DEBUG

#include "io.h"

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <unistd.h>
#include <math.h>

// not the most efficient way...
// but works for now
#define MIN(x,y) (x < y ? x : y)

FILE* x_file;
FILE* y_file;
FILE* b_file;

double* x[2];
double* y[2];
double* b[2];


double *x_compute, *y_compute, *b_compute;

sem_t sem_io;
sem_t sem_comp;

typedef struct problem_args_t {
  int p;
  int n;
  int m;
  int t;
  int x_b;
  int y_b;
  int m_indexed;
  int t_indexed;
} problem_args;

problem_args in;

void write_test_matrices(problem_args *args) {
  double* out;
  int i, j;
  int len_x = args->p*args->n*args->m;
  int len_y = args->n*args->t;
  printf("creating test matrices\n");
  out = (double*)malloc(len_x*sizeof(double));
  for (j = 0; j < args->m; j++) {
    for (i = 0; i < args->p*args->n; i++) {
      out[i + args->p*args->n*j] = j+1;
    }
  }
  write_double(out, x_file, args->p*args->n, args->m, 0);
  free(out);

  out = (double*)malloc(len_y*sizeof(double));
  for (j = 0; j < args->t; j++) {
    for (i = 0; i < args->n; i++) {
      out[i + j*args->n] = j+1;
    }
  }
  write_double(out, y_file, args->n, args->t, 0);
  free(out);
}

int return_buffer_index(double** buffers, int size, double* cur) {
  int i;
  for(i = 0; i < size; i++) {
    if (buffers[i] == cur) {
      return i;
    }
  }
  return -1;
}

void swap_buffers(double** b1, double** b2) {
  double* temp;
  temp = *b1;
  *b1 = *b2;
  *b2 = temp;
}

void read_x(double* buf, int index, const problem_args* args) {
#ifdef DEBUG
  printf("read_x:\n\tx[%d](x=%d)\n", return_buffer_index(x, 2, buf), index);  
#endif
  read_double(buf, x_file, args->p*args->n,
              MIN(args->x_b, args->m - args->x_b*index),
              index);
}

void read_y(double* buf, int index, const problem_args* args) {
#ifdef DEBUG
  printf("read_y:\n\ty[%d](y=%d)\n", return_buffer_index(y, 2, buf), index);  
#endif
  read_double(buf, y_file, args->n, 
              MIN(args->y_b, args->t - args->y_b*index), 
              index);
}

void write_b(double* buf, int s, int r, const problem_args* args) {
  int j, y_inc, x_inc;
  y_inc = MIN(args->y_b, args->t - args->y_b*r);
  x_inc = MIN(args->x_b, args->m - args->x_b*s);
#ifdef DEBUG
  printf("write:\n\tb[%d](x=%d)(y=%d)\n", return_buffer_index(b, 2, buf), s, r); 
  printf("\tstarting at %d\n", s*x_inc*args->p+r*args->p*args->m);
#endif
  write_double(buf, b_file, args->p, 
	       x_inc, s*x_inc*args->p+r*args->p*args->m);
}

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

  read_x(x_cur, 0, args);
  read_y(y_cur, 0, args);
  sem_post(&sem_comp);

  for (r = 0; r < i; r++) {
    swap_buffers(&y_cur, &y_next);
    read_y(y_cur, (r+1) % i, args);
    for (s = 0; s < j; s++) {
      swap_buffers(&x_cur, &x_next);
      read_x(x_cur, (s+1) % j, args);
      sem_post(&sem_comp);
      sem_wait(&sem_io);
      swap_buffers(&b_prev, &b_cur);
      write_b(b_prev, s, r, args);
    }
  }
  pthread_exit(NULL);
}

#define ITEM(arr, x, x_s, y, y_s, z) arr[x + x_s*y + x_s*y_s*z]
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
      sem_wait(&sem_comp);
#ifdef DEBUG
      printf("compute(x_cur)\n:");
      print_buffer(x_cur, args->p*args->n*args->x_b);
      printf("compute(y_cur):\n");
      print_buffer(y_cur, args->n*args->y_b);
#endif
      computation(x_cur, y_cur, b_cur, args);
      printf("compute:\n\tx[%d](x=%d)\n\ty[%d](y=%d)\n\tb[%d]\n", 
	     return_buffer_index(x, 2, x_cur), s,
	     return_buffer_index(y, 2, y_cur), r,
	     return_buffer_index(b, 2, b_cur));
      printf("compute(b_cur):\n");
      print_buffer(b_cur, args->p*args->y_b*args->x_b);
      swap_buffers(&x_cur, &x_next);
      swap_buffers(&b_prev, &b_cur);
      sem_post(&sem_io);
    }
    swap_buffers(&y_cur, &y_next);
  }
  pthread_exit(NULL);

}

void print_output(problem_args *args) {
  double *buf;
  int i, j;
  buf = (double*) malloc(args->p*args->m*args->t*sizeof(double));
  read_double( buf, b_file, args->p*args->m*args->t, 1, 0);
  printf("printing output:\n");
  print_buffer( buf, args->p*args->m*args->t);
}
int main(int argc, char* argv[]) {
  int rc;
  int i;
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
  
  write_test_matrices(&in);
  
  x[0] =  (double*)malloc(in.p * in.n * in.x_b * sizeof(double));
  x[1] = (double*)malloc(in.p * in.n * in.x_b * sizeof(double));
  y[0] =  (double*)malloc(in.n * in.y_b * sizeof(double));
  y[1] = (double*)malloc(in.n * in.y_b * sizeof(double));
  b[0] = (double*)malloc(in.p * in.x_b * in.y_b *sizeof(double));
  b[1] =  (double*)malloc(in.p * in.x_b * in.y_b *sizeof(double));
  
  sem_init(&sem_io, 0, 0);
  sem_init(&sem_comp, 0, 0);

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

  fflush(b_file);

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
