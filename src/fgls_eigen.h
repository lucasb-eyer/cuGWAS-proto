#ifndef FGLS_EIGEN_H
#define FGLS_EIGEN_H

#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
#define STR_BUFFER_SIZE 256

#include "timing.h"

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
#endif // TIMING

} FGLS_eigen_t;

FGLS_eigen_t FGLS_eigen_config;

int  fgls_eigen( FGLS_eigen_t *cf );
int  preloop(double *Phi, double *Z, double *W);
void swap_buffers(double** b1, double** b2);

#endif // FGLS_H
