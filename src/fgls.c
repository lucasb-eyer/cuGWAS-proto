#include "fgls.h"
#include "io.h"

#include <stdio.h>
#include <stdlib.h>

void print_output(const problem_args *args) {
  double *buf;
  if(!b_file) {
    printf("b_file not initialized. Exiting...\n");
    exit(-1);
  }
  buf = (double*) malloc(args->p*args->m*args->t*sizeof(double));
  read_double(buf, b_file, args->p*args->m*args->t, 1, 0);
  printf("printing output:\n");
  print_buffer( buf, args->p*args->m*args->t);
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
  read_double(buf, x_file, args->p*args->n,
              MIN(args->x_b, args->m - args->x_b*index),
              index);
}

void read_y(double* buf, int index, const problem_args* args) {
  if(!y_file) {
    printf("y_file not initialized. Exiting...\n");
    exit(-1);
  }
  read_double(buf, y_file, args->n, 
              MIN(args->y_b, args->t - args->y_b*index), 
              index);
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


void write_b(double* buf, int s, int r, const problem_args* args) {
  if(!b_file) {
    printf("b_file not initialized. Exiting...\n");
    exit(-1);
  }
  int y_inc, x_inc;
  y_inc = MIN(args->y_b, args->t - args->y_b*r);
  x_inc = MIN(args->x_b, args->m - args->x_b*s);
  write_double(buf, b_file, args->p, 
	       x_inc, s*x_inc*args->p+r*args->p*args->m);
}


void write_test_matrices(FILE* x_file, FILE* y_file, problem_args *args) {
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
