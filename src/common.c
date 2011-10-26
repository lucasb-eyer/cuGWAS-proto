#include "common.h"


void initialize_config(
		FGLS_config_t *cf,
		int n, int p, int m, int t, int wXL, int wXR,
        int x_b, int y_b, int num_threads,
		char *Phi_path, char *h2_path, char *sigma2_path,
		char *XL_path, char *XR_path, char *Y_path,
		char *B_path, char *V_path
)
{
	// Problem dimensions
	cf->n = n;
	cf->p = p;
	cf->m = m;
	cf->t = t;
	cf->wXR = wXR;
	cf->wXL = wXL;
	// Algorithm parameters
	cf->x_b = x_b;
	cf->y_b = y_b;
	cf->NUM_COMPUTE_THREADS = num_threads;
	// In/Out Files
	strncpy( cf->Phi_path,    Phi_path,    STR_BUFFER_SIZE );
	strncpy( cf->h_path,      h2_path,     STR_BUFFER_SIZE );
	strncpy( cf->sigma_path , sigma2_path, STR_BUFFER_SIZE );
	strncpy( cf->XL_path,     XL_path,     STR_BUFFER_SIZE );
	strncpy( cf->XR_path,     XR_path,     STR_BUFFER_SIZE );
	strncpy( cf->Y_path,      Y_path,      STR_BUFFER_SIZE );
	strncpy( cf->B_path,      B_path,      STR_BUFFER_SIZE );
	strncpy( cf->V_path,      V_path,      STR_BUFFER_SIZE );
	// Temporary files
	snprintf( cf->ZtXL_path, STR_BUFFER_SIZE, "%s.tmp", XL_path );
	snprintf( cf->ZtXR_path, STR_BUFFER_SIZE, "%s.tmp", XR_path );
	snprintf( cf->ZtY_path,  STR_BUFFER_SIZE, "%s.tmp", Y_path );

	return;
}

void swap_buffers(double** b1, double** b2) 
{
	double* tmp;

	tmp = *b1;
	*b1 = *b2;
	*b2 = tmp;
}

void swap_aiocb(struct aiocb ***x, struct aiocb ***y)
{
	struct aiocb **tmp = *x;
	*x = *y;
	*y = tmp;
}

