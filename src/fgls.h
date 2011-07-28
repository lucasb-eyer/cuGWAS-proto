#ifndef FGLS_H
#define FGLS_H

#include <stdio.h>
#include <sys/time.h>
#include <time.h>

// not the most efficient way...
// but works for now
#define MIN(x,y) (x < y ? x : y)

// get item at (x,y,z) of array 'arr'
#define ITEM(arr, x, x_s, y, y_s, z) arr[x + x_s*y + x_s*y_s*z]

#define DEBUG 0
#define TIMING 1

FILE* x_file;
FILE* y_file;
FILE* h_file;
FILE* phi_file;
FILE* b_file;

typedef struct timing_t {
  long compute_time;
  long io_time;
  long comp_mutex_wait_time;
  long io_mutex_wait_time;
} timing;

typedef struct problem_args_t {
  int p;
  int n;
  int m;
  int t;
  int x_b;
  int y_b;
  int m_indexed;
  int t_indexed;
  //  double h;
  timing* time;
} problem_args;

double compare(double *a, double *b, int size);
long get_diff_ms(struct timeval *s, struct timeval *e);
void print_output(FILE* f, const problem_args* args);
void swap_buffers(double** b1, double** b2);
void read_x(double* buf, int index, const problem_args* args);
void read_phi(double* buf, int index, const problem_args* args);
void read_y(double* buf, int index, const problem_args* args);
void read_h(double* buf, int index, const problem_args* args);
int  return_buffer_index(double** buffers, int size, double* cur);
void write_b(double* buf, int s, int r, const problem_args* args);
void write_x(double* buf, int s, const problem_args* args);
void write_y(double* buf, int r, const problem_args* args);
void write_test_matrices(FILE* x_file, FILE* y_file, problem_args *args);

#endif // FGLS_H
