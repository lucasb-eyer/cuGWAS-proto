#ifndef FGLS_EIGEN_H
#define FGLS_EIGEN_H

#include "common.h"

typedef struct 
{
	// In-core operands
	double *Z;
	double *W;
	double *h2;
	double *sigma2;
	double *XL_orig;
	double *ZtXL;
	// Out-of-core operands
	FILE *fp_Y;
	FILE *fp_ZtY;
	// Output
	double *res_sigma;
} ooc_res_sigma_t;

typedef struct 
{
	// Dimensions
	int m, n, k;
	long int n_cols_per_buff;
	// Operands
	double *Z;
	FILE *fp_in;
	FILE *fp_out;
	double *in[2];
	double *out[2];
} ooc_gemm_t;

typedef struct {
	// In-core operands
	/*double *Z;*/
	double *W;
	double *h;
	double *sigma;
	double *res_sigma;
	double *alpha;
	double *beta;
	double *Winv;

	double *XL[2];
	double *XL_orig;
	double *B_t;
	double *V_tl;

	// Out-of-core operands
	double *X[NUM_BUFFERS_PER_THREAD];
	double *Y[NUM_BUFFERS_PER_THREAD];
	double *B[NUM_BUFFERS_PER_THREAD];
	double *V[NUM_BUFFERS_PER_THREAD];

	FILE *XR_fp;
	FILE *Y_fp;
	FILE *B_fp;
	FILE *V_fp;

	// Configuration
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
