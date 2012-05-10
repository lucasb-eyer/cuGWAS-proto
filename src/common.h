#ifndef FGLS_COMMON_H
#define FGLS_COMMON_H

#define DEBUG 1
#define TIMING 0
#define VAMPIR 0

#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
#define STR_BUFFER_SIZE 256
#define NUM_BUFFERS_PER_THREAD 2

#include "aio.h"

typedef struct
{
  char Phi_path[STR_BUFFER_SIZE];
  char h_path[STR_BUFFER_SIZE];
  char sigma_path[STR_BUFFER_SIZE];
  char XL_path[STR_BUFFER_SIZE];
  char XR_path[STR_BUFFER_SIZE];
  char ZtXL_path[STR_BUFFER_SIZE];
  char ZtXR_path[STR_BUFFER_SIZE];
  char Y_path[STR_BUFFER_SIZE];
  char ZtY_path[STR_BUFFER_SIZE];
  char B_path[STR_BUFFER_SIZE];
  char V_path[STR_BUFFER_SIZE];

  int n;
  int p; // total width of X
  int m;
  int t;
  int wXL; // width of XL
  int wXR; // width of XR
  int x_b;
  int y_b;

  int NUM_COMPUTE_THREADS;
} FGLS_config_t;

void initialize_config(
        FGLS_config_t *cf,
        int n, int p, int m, int t, int wXL, int wXR,
        int x_b, int y_b, int num_threads,
        char *Phi_path, char *h2_path, char *sigma2_path,
        char *XL_path, char *XR_path, char *Y_path,
        char *B_path, char *V_path
);

void swap_buffers(double **b1, double **b2);
void swap_aiocb(const struct aiocb ***x, const struct aiocb ***y);

int read_clock(struct timeval *t);
int elapsed_time(struct timeval *start, struct timeval *end);

#define fgls_malloc(size) fgls_malloc_impl(__FILE__, __LINE__, size)
void * fgls_malloc_impl( const char* file, long line, size_t size );
void error_msg(char *msg, int abort);

#endif // FGLS_COMMON_H
