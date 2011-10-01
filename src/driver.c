#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "io.h"
#include "timing.h"
#include "fgls_eigen.h"

double compare(double *a, double *b, int size);

int main( int argc, char *argv[] ) 
{
  char *dir;
  FGLS_eigen_t *cf = &FGLS_eigen_config;

  if (argc != 2) 
  {
    fprintf(stderr, "Usage: %s <data_set_directory>\n", argv[0]);
	fprintf(stderr, "\tnote: data_set_directory must contain files:\n");
	fprintf(stderr, "X.in\nY.in\nPhi.in\nH.in\nSig.in\n(B_exp.out if checking answer)\n");
    exit(EXIT_FAILURE);
  }

  dir = argv[1];
  printf("Please enter parameters\n");
  printf("n: ");                 scanf("%d", &cf->n);
  printf("p: ");                 scanf("%d", &cf->p);
  printf("m: ");                 scanf("%d", &cf->m);
  printf("t: ");                 scanf("%d", &cf->t);
  printf("width of XL: ");       scanf("%d", &cf->wXL);
  printf("x_b: ");               scanf("%d", &cf->x_b);
  printf("# compute threads: "); scanf("%d", &cf->NUM_COMPUTE_THREADS);
  cf->wXR = cf->p - cf->wXL;

#if TIMING
  cf->time = (timing*)malloc(sizeof(timing));
  struct timeval start, end;
  gettimeofday(&start, NULL);
#endif // TIMING

  sprintf(cf->XR_path,   "%s/XR.in", dir);
  sprintf(cf->XL_path,   "%s/XL.in", dir);
  sprintf(cf->ZtXL_path, "%s/XL.tmp", dir);
  sprintf(cf->ZtXR_path, "%s/XR.tmp", dir);
  sprintf(cf->Y_path,    "%s/Y.in", dir);
  sprintf(cf->ZtY_path,  "%s/Y.tmp", dir);
  sprintf(cf->Phi_path,  "%s/Phi.in", dir);
  sprintf(cf->h_path,    "%s/H.in", dir);
  sprintf(cf->sigma_path,"%s/Sig.in", dir);
  sprintf(cf->B_path,    "%s/B.out", dir);

  fgls_eigen( cf );

#if TIMING
  gettimeofday(&end, NULL);
  long total = get_diff_ms(&start, &end);;
#endif // TIMING
#if DEBUG
  char str_buf[STR_BUFFER_SIZE];
  double *b_mine, *b_exp;
  FILE *b_mine_f, *b_exp_f;
  b_mine = (double*)malloc(cf->m*cf->t*cf->p*sizeof(double));
  b_exp = (double*)malloc(cf->m*cf->t*cf->p*sizeof(double));

  sprintf(str_buf, "%s/B.out", dir);
  b_mine_f = fopen(str_buf, "rb");
  if(!b_mine_f) {
    printf("ERROR opening %s\n", str_buf);
    return -1;
  }

  sprintf(str_buf, "%s/B_exp.out", dir);
  b_exp_f = fopen(str_buf, "rb");
  if(!b_exp_f) {
    printf("ERROR opening %s\n", str_buf);
    return -1;
  }
  double max;
  sync_read(b_mine, b_mine_f, cf->p*cf->m*cf->t, 0);
  sync_read(b_exp,  b_exp_f,  cf->p*cf->m*cf->t, 0);
  /*printf("out[0]: %12e\n", b_mine[0]);*/
  /*printf("exp[0]: %12e\n", b_exp[0]);*/
  /*printf("out[1]: %12e\n", b_mine[1]);*/
  /*printf("exp[1]: %12e\n", b_exp[1]);*/
  /*printf("out[3000]: %12e\n", b_mine[3000]);*/
  /*printf("exp[3000]: %12e\n", b_exp[3000]);*/
  /*printf("out[6000]: %12e\n", b_mine[6000]);*/
  /*printf("exp[6000]: %12e\n", b_exp[6000]);*/
  /*printf("out[9000]: %12e\n", b_mine[9000]);*/
  /*printf("exp[9000]: %12e\n", b_exp[9000]);*/
  max = compare(b_mine, b_exp, cf->m*cf->p*cf->t);
  printf("Max elemental diff: %lf\n", max);
#endif // DEBUG

#if TIMING
  if (!cf->time) {
    return 0;
  }
  printf("-- timing_info [    x_b    total(ms) ] = [      IO    |    COMP   |  IO_wait  | COMP_wait ];\n");
  printf("-- timing_info [ %6d    %6ld    ] = [ %10ld |%10ld |%10ld |%10ld ];\n", 
		  cf->x_b, 
		  total, 
		  cf->time->io_time, 
		  cf->time->compute_time, 
          cf->time->io_mutex_wait_time, 
		  cf->time->comp_mutex_wait_time);
#endif // TIMING
  return 0;
}

double compare(double* a, double* b, int size) {
  double out = 0.0;
  int i;
  for (i = 0; i < size; i++) {
    if((i % 1000) == 0 && fabs(a[i] - b[i]) > 1e-13)
		printf("Difference at %d: %e [%f - %f]\n", i, fabs(a[i] - b[i]), a[i], b[i]);
    if(out < fabs(a[i] - b[i]))
      out = fabs(a[i] - b[i]);
  }
  return out;
}

