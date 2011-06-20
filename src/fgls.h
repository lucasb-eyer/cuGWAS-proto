#ifndef FGLS_H
#define FGLS_H

#include <stdio.h>

// not the most efficient way...
// but works for now
#define MIN(x,y) (x < y ? x : y)

// get item at (x,y,z) of array 'arr'
#define ITEM(arr, x, x_s, y, y_s, z) arr[x + x_s*y + x_s*y_s*z]

// if you are working on OSX...
#define USE_NAMED_SEMAPHORES 0
#define DEBUG 1

FILE* x_file;
FILE* y_file;
FILE* b_file;

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

void print_output(const problem_args* args);
void swap_buffers(double** b1, double** b2);
void read_x(double* buf, int index, const problem_args* args);
void read_y(double* buf, int index, const problem_args* args);
int  return_buffer_index(double** buffers, int size, double* cur);
void write_b(double* buf, int s, int r, const problem_args* args);
void write_test_matrices(FILE* x_file, FILE* y_file, problem_args *args);

#endif // FGLS_H
