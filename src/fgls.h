#ifndef FGLS_H
#define FGLS_H

#define DEBUG 0
#define TIMING 0

#include <stdio.h>

#if TIMING
#include <sys/time.h>
#include <time.h>
#endif // TIMING

#if TIMING
#define BEGIN_TIMING() \
    gettimeofday(&start, NULL);
#define END_TIMING(var) \
  gettimeofday(&end, NULL); \
  var += get_diff_ms(&start, &end);
#define DEF_TIMING() \
  struct timeval start, end;
#else // TIMING
#define BEGIN_TIMING()
#define END_TIMING(var)
#define DEF_TIMING()
#endif // TIMING

// not the most efficient way...
// but works for now
#define MIN(x,y) (x < y ? x : y)

// get item at (x,y,z) of array 'arr'
#define ITEM(arr, x, x_s, y, y_s, z) arr[x + x_s*y + x_s*y_s*z]

#define STR_BUFFER_SIZE 256



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
#if TIMING
  timing* time;
#endif TIMING
} problem_args;

#if TIMING
long get_diff_ms(struct timeval *s, struct timeval *e);
#endif // TIMING

double compare(double *a, double *b, int size);
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
