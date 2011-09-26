#ifndef FGLS_H
#define FGLS_H

#define DEBUG 0
#define TIMING 0
#define VAMPIR 1

#include <stdio.h>

/*#if TIMING*/
#include <sys/time.h>
#include <time.h>
/*#endif // TIMING*/

#if TIMING
#define BEGIN_TIMING() \
	gettimeofday(&start, NULL);
#define END_TIMING(var)     \
	gettimeofday(&end, NULL);   \
var += get_diff_ms(&start, &end);
#define END_TIMING2(var)     \
	gettimeofday(&end, NULL);   \
/*printf("Time: %ld\n", get_diff_ms(&start, &end)); \*/ \
var += get_diff_us(&start, &end);
#define DEF_TIMING() \
	struct timeval start, end;
#else // TIMING
#define BEGIN_TIMING()
#define END_TIMING(var)
#define END_TIMING2(var)
#define DEF_TIMING()
#endif // TIMING

/*#if TIMING*/
/*#define BEGIN_TIMING() \*/
/*clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);*/
/*#define END_TIMING(var)     \*/
/*clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end); \*/
/*var += get_diff_ns(&start, &end);*/
/*#define END_TIMING2(var)     \*/
/*gettimeofday(&end, NULL);   \*/
/*printf("Time: %ld\n", get_diff_ms(&start, &end)); \*/
/*var += get_diff_us(&start, &end);*/
/*#define DEF_TIMING() \*/
/*struct timespec start, end;*/
/*#else // TIMING*/
/*#define BEGIN_TIMING()*/
/*#define END_TIMING(var)*/
/*#define DEF_TIMING()*/
/*#endif // TIMING*/

// not the most efficient way...
// but works for now
#define MIN(x,y) (x < y ? x : y)

// get item at (x,y,z) of array 'arr'
/*#define ITEM(arr, x, x_s, y, y_s, z) arr[x + x_s*y + x_s*y_s*z]*/

#define STR_BUFFER_SIZE 256

typedef struct timing_t {
  long compute_time;
  long io_time;
  long comp_mutex_wait_time;
  long io_mutex_wait_time;
} timing;

typedef struct
{
  char XL_path[STR_BUFFER_SIZE];
  char XR_path[STR_BUFFER_SIZE];
  char ZtXL_path[STR_BUFFER_SIZE];
  char ZtXR_path[STR_BUFFER_SIZE];
  char Y_path[STR_BUFFER_SIZE];
  char ZtY_path[STR_BUFFER_SIZE];
  char Phi_path[STR_BUFFER_SIZE];
  char h_path[STR_BUFFER_SIZE];
  char sigma_path[STR_BUFFER_SIZE];
  char B_path[STR_BUFFER_SIZE];

  int n;
  int p; // total width of X
  int m;
  int t;
  int wXL; // width of XL
  int wXR; // width of XR
  int x_b;
  /*int y_b;*/

  int NUM_COMPUTE_THREADS;
  /*int NUM_BUFFERS_PER_THREAD;*/

#if TIMING
  timing* time;
#endif TIMING

} FGLS_eigen_t;

FGLS_eigen_t FGLS_eigen_config;

/*typedef struct problem_args_t {*/
/*int id;*/
/*} problem_args;*/

long get_diff_ms(struct timeval *s, struct timeval *e);
#if TIMING
long get_diff_us(struct timeval *s, struct timeval *e);
long get_diff_ns(struct timespec *s, struct timespec *e);
#endif // TIMING

double compare(double *a, double *b, int size);
void   swap_buffers(double** b1, double** b2);
/*void   read_x(double* buf, int index, const problem_args* args);*/
/*void   read_phi(double* buf, int index, const problem_args* args);*/
/*void   read_y(double* buf, int index, const problem_args* args);*/
/*void   read_h(double* buf, int index, const problem_args* args);*/
/*int    return_buffer_index(double** buffers, int size, double* cur);*/
/*void   write_b(double* buf, int s, int r, const problem_args* args);*/
/*void   write_x(double* buf, int s, const problem_args* args);*/
/*void   write_y(double* buf, int r, const problem_args* args);*/

#endif // FGLS_H
