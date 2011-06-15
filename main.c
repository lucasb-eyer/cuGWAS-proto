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

double* x[2];
double* y[2];
double* b[2];

double *x_cur;
double *x_next;
double *y_cur;
double *y_next;
double *b_prev;
double *b_cur;
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
  int i;
  int len_x = args->p*args->n*args->m;
  int len_y = args->n*args->t;
  printf("creating test matrices\n");
  out = (double*)malloc(len_x*sizeof(double));
  for( i = 0; i < len_x; i++) {
    out[i] = i;
  }
  write_double(out, "x.in", args->p*args->n, args->m, 0);
  free(out);

  out = (double*)malloc(len_y*sizeof(double));
  for( i = 0; i < len_y; i++) {
    out[i] = i;
  }
  write_double(out, "y.in", args->n, args->t, 0);
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

void read_x(int index, const problem_args* args) {
  printf("read_x:\n\tx[%d](x=%d)\n", return_buffer_index(x, 2, x_next), index);  
  read_double(x_next, "x.in", args->p*args->n,
              MIN(args->x_b, args->m - args->x_b*index),
              index);
}

void read_y(int index, const problem_args* args) {
  printf("read_y:\n\ty[%d](y=%d)\n", return_buffer_index(y, 2, y_next), index);  
  read_double(y_next, "y.in", args->n, 
              MIN(args->y_b, args->t - args->y_b*index), 
              index);
}

void write_b(int s, int r, const problem_args* args) {
  int z, y_inc, x_inc;
  y_inc =  MIN(args->y_b, args->t - args->y_b*r);
  x_inc = MIN(args->x_b, args->m - args->x_b*s);
  printf("write:\n\tb[%d](x=%d)(y=%d)\n", return_buffer_index(b, 2, b_prev), s, r);  
  for( z = r; z < args->t; z += y_inc) {
    write_double(b_prev, "b.out", args->p, 
		 x_inc, s);
  }
}

void* io(void* in) {
  problem_args* args = (problem_args*)in;
  int s, r;
  int i = args->t_indexed;
  int j = args->m_indexed;

  // weird starting condition
  // for now, read into both buffers
  // need to fix in the future
  read_x(0, args);
  read_y(0, args);
  swap_buffers(&y_cur, &y_next);
  swap_buffers(&x_cur, &x_next);
  read_x(0, args);
  read_y(0, args);
  sem_post(&sem_comp);

  for (r = 0; r < i; r++) {
    swap_buffers(&y_cur, &y_next);
    read_y((r+1) % i, args);
    for (s = 0; s < j; s++) {
      swap_buffers(&x_cur, &x_next);
      read_x((s+1) % j, args);
      sem_post(&sem_comp);
      sem_wait(&sem_io);
      swap_buffers(&b_prev, &b_cur);
      write_b(s, r, args);
    }
  }
  pthread_exit(NULL);
}

#define ITEM(arr, x, x_s, y, y_s, z) arr[y_s*x + y + x_s*y_s*z]
// compute a (dot) b = c 
//    a - args.p(h) x args.n(k) x args.x_b(i)
//    b - args.n(k) x 1 x args.y_b(j)
//    c - args.p(h) x args.x_b(i) x args.y_b(j)
void computation( double* a, double* b, double* c, problem_args* args) {
  int i, j, h, k;
  double sum;
  printf("inside computation:\n");
  print_buffer( a, 16);
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

  int i = args->t_indexed;
  int j = args->m_indexed;

  for (r = 0; r < i; r++) {
    for (s = 0; s < j; s ++) {
      sem_wait(&sem_comp);
      x_compute = x_cur;
      y_compute = y_cur;
      b_compute = b_cur;
      sem_post(&sem_io);
      computation(x_compute, y_compute, b_compute, args);
      printf("compute:\n\tx[%d](x=%d)\n\ty[%d](y=%d)\n\tb[%d]\n", 
	     return_buffer_index(x,2,x_compute), s,
	     return_buffer_index(y,2,y_compute), r,
	     return_buffer_index(b,2,b_compute));
      sleep(1);
    }
  }
  pthread_exit(NULL);

}

void print_output(problem_args *args) {
  double *buf;
  buf = (double*) malloc(args->p*args->m*args->t*sizeof(double));
  read_double( buf, "b.out", args->p*args->m*args->t, 1, 0);
  printf("printing output:\n");
  print_buffer( buf, args->p*args->m*args->t);
}
int main() {
  int rc;
  int i;
  pthread_t io_thread;
  pthread_t compute_thread;

  in.m = 1;
  in.t = 1;

  in.x_b = 1;
  in.y_b = 1;
  in.n = 4;
  in.p = 4;
  in.m_indexed = (int) ((double) in.m/in.x_b);
  in.t_indexed = (int) ((double) in.t/in.y_b);
  
  write_test_matrices(&in);
  
  x_cur =  (double*)malloc(in.p * in.n * in.x_b * sizeof(double));
  x_next = (double*)malloc(in.p * in.n * in.x_b * sizeof(double));
  y_cur =  (double*)malloc(in.n * in.y_b * sizeof(double));
  y_next = (double*)malloc(in.n * in.y_b * sizeof(double));
  b_prev = (double*)malloc(in.p * in.x_b * in.y_b *sizeof(double));
  b_cur =  (double*)malloc(in.p * in.x_b * in.y_b *sizeof(double));
  
  x[0] = x_cur;
  x[1] = x_next;
  y[0] = y_cur;
  y[1] = y_next;
  b[0] = b_prev;
  b[1] = b_cur;
  printf("0x%x\n", sizeof(double)*in.p*in.x_b);
  for(i = 0; i < 2; i++) {
    printf("x[%d] = 0x%x\n", i, x[i]);
  }
  for(i = 0; i < 2; i++) {
    printf("y[%d] = 0x%x\n", i, y[i]);
  }
  for(i = 0; i < 2; i++) {
    printf("b[%d] = 0x%x\n", i, b[i]);
  }
  sem_init(&sem_io, 0, 0);
  sem_init(&sem_comp, 0, 0);

  rc = pthread_create( &io_thread, NULL, io, (void*)&in);
  if (rc){
    printf("error: return code from pthread_create() is %d\n", rc);
    pthread_exit(NULL);
  }
  rc = pthread_create( &compute_thread, NULL, compute, (void*)&in);
  if (rc){
    printf("error: return code from pthread_create() is %d\n", rc);
    pthread_exit(NULL);
  }

  void* retval;

  pthread_join(io_thread, &retval);
  pthread_join(compute_thread, &retval);

  print_output(&in);

  pthread_exit(NULL);
  return 0;
}
