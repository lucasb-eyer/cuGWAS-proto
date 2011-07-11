#include "fgls.h"
#include "io.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

double compare(double* a, double* b, int size) {
  double out = 0.0;
  int i;
  printf("a:\t\tb:\n");
  for (i = 0; i < size; i++) {
    printf("%lf\t%lf\n", a[i], b[i]);
    if(out < fabs( a[i] - b[i]))
      out = fabs(a[i] - b[i]);
  }
  return out;
}

long get_diff_ms(struct timeval *start_time, struct timeval *end_time) {
  long seconds = end_time->tv_sec - start_time->tv_sec;
  long useconds = end_time->tv_usec - start_time->tv_usec;
  return ((seconds) * 1000 + useconds/1000.0) +0.5;
}
void print_output(FILE* b_file, const problem_args *args) {
  double *buf;
  if(!b_file) {
    printf("b_file not initialized. Exiting...\n");
    exit(-1);
  }
  buf = (double*) malloc(args->p*args->m*args->t*sizeof(double));
  read(buf, b_file, args->p*args->m*args->t, 0);
  printf("printing output:\n");
  print_buffer( buf, args->p*args->m*args->t);
  free(buf);
}

void swap_buffers(double** b1, double** b2) {
  double* temp;
  temp = *b1;
  *b1 = *b2;
  *b2 = temp;
}

void read_x(double* buf, int index, const problem_args* args) {
  if(!x_file) {
    printf("x_file not initialized. Exiting...\n");
    exit(-1);
  }
  int x_inc = MIN(args->x_b, args->m - args->x_b*index);
  read(buf, x_file, args->p*args->n*x_inc, index*args->x_b*args->p*args->n);
}

void read_phi(double* buf, int index, const problem_args* args) {
  if(!phi_file) {
    printf("phi_file not initialized. Exiting...\n");
    exit(-1);
  }
  read(buf, phi_file, args->n*args->n, 0);
}

void read_y(double* buf, int index, const problem_args* args) {
  if(!y_file) {
    printf("y_file not initialized. Exiting...\n");
    exit(-1);
  }
  int y_inc = MIN(args->y_b, args->t - args->y_b*index); 
  read(buf, y_file, args->n*y_inc, index*args->y_b*args->n);
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


void write_b(double* buf, int r, int s, const problem_args* args) {
  if(!b_file) {
    printf("b_file not initialized. Exiting...\n");
    exit(-1);
  }
  //  printf("r: %d s: %d\n", r, s );
  int y_inc, x_inc, j, buffer_index, file_index;
  y_inc = MIN(args->y_b, args->t - args->y_b*s);
  x_inc = MIN(args->x_b, args->m - args->x_b*r);
  for (j = 0; j < y_inc; j++) {
    buffer_index = args->x_b*args->p*j;
    file_index = args->m*args->p*(args->y_b*s + j) + r*args->x_b*args->p;
    write(&buf[buffer_index], b_file, args->p*x_inc, file_index);
  }
  
}

/*
void write_test_matrices(FILE* x_file, FILE* y_file, problem_args *args) {
  double* out;
  int i, j;
  int len_x = args->p*args->n*args->m;
  int len_y = args->n*args->t;
#ifdef DEBUG
  printf("creating test matrices\n");
#endif
  out = (double*)malloc(len_x*sizeof(double));
  for (j = 0; j < args->m; j++) {
    for (i = 0; i < args->p*args->n; i++) {
      out[i + args->p*args->n*j] = j+1;
    }
  }
  write(out, x_file, args->p*args->n, args->m, 0);
  free(out);

  out = (double*)malloc(len_y*sizeof(double));
  for (j = 0; j < args->t; j++) {
    for (i = 0; i < args->n; i++) {
      out[i + j*args->n] = j+1;
    }
  }
  write(out, y_file, args->n, args->t, 0);
  free(out);
}
*/
