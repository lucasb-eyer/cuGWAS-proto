#ifndef FGLS_EIGEN_H
#define FGLS_EIGEN_H

#include <semaphore.h>
#include "common.h"

#define NUM_BUFFERS_PER_THREAD 2

typedef struct 
{
	double *Z;
	FILE *fp_in;
	FILE *fp_out;
	int m, n, k;
	long int n_cols_per_buff;
	double *in[2];
	double *out[2];
	sem_t sem_io;
	sem_t sem_comp;
} ooc_gemm_t;

typedef struct {
	FILE *XR_fp;
	FILE *Y_fp;
	FILE *B_fp;
	FILE *V_fp;
	double *XL[2];
	double *XL_b;
	double *XLtXL;

	double *X[NUM_BUFFERS_PER_THREAD];
	double *Y[NUM_BUFFERS_PER_THREAD];
	double *B[NUM_BUFFERS_PER_THREAD];
	double *V[NUM_BUFFERS_PER_THREAD];

	double *h;
	double *sigma;
	double *W;
	double *alpha;
	double *beta;
	double *Winv;
	double *xtSx;

	FGLS_config_t *cf;

	int id;
} ooc_loops_t;

int fgls_eigen(
		int n, int p, int m, int t, int wXL, int wXR,
        int x_b, int y_b, int num_threads,
		char *Phi_path, char *h2_path, char *sigma2_path,
		char *XL_path, char *XR_path, char *Y_path,
		char *B_path, char *V_path
);

#endif // FGLS_EIGEN_H
