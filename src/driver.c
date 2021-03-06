#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "io.h"
#include "timing.h"
#include "fgls_chol.h"
#ifdef FGLS_WITH_GPU
    #include "fgls_chol_gpu.h"
#endif
#include "fgls_eigen.h"

double compare(double *a, double *b, int size);

int main( int argc, char *argv[] ) 
{
  char *dir, var;
  int i, nrep;
  FGLS_config_t cf;
  struct timeval start, end;

  if (argc != 4) 
  {
    fprintf(stderr, "Usage: %s <variant - [e|m]> <data_set_directory> <#repetitions>\n", argv[0]);
    fprintf(stderr, "\tnote: data_set_directory must contain files:\n");
    fprintf(stderr, "- Phi.in\n- H.in\n- Sig.in\n- XL.in\n- XR.in\n- Y.in\n");
    exit(EXIT_FAILURE);
  }

  var  = argv[1][0];
  dir  = argv[2];
  nrep = atoi(argv[3]);
  /*printf("Please enter parameters\n");*/
  /*printf("n: ");*/
  scanf("%d", &cf.n);
  /*printf("p: ");*/
  scanf("%d", &cf.p);
  /*printf("m: ");   */
  scanf("%d", &cf.m);
  /*printf("t: "); */
  scanf("%d", &cf.t);
  /*printf("width of XL: ");*/
  scanf("%d", &cf.wXL);
  /*printf("x_b: "); */
  scanf("%d", &cf.x_b);
  /*printf("y_b: "); */
  scanf("%d", &cf.y_b);
  /*printf("# compute threads: ");*/
  scanf("%d", &cf.NUM_COMPUTE_THREADS);
  cf.wXR = cf.p - cf.wXL;

  gettimeofday(&start, NULL);

  for ( i = 0; i < nrep; i++ )
  {
      sprintf(cf.XL_path,   "%s/XL.in", dir);
      sprintf(cf.XR_path,   "%s/XR.in", dir);
      /*sprintf(cf.ZtXL_path, "%s/XL.tmp", dir);*/
      /*sprintf(cf.ZtXR_path, "%s/XR.tmp", dir);*/
      sprintf(cf.Y_path,    "%s/Y.in", dir);
      sprintf(cf.ZtY_path,  "%s/Y.tmp", dir);
      sprintf(cf.Phi_path,  "%s/Phi.in", dir);
      sprintf(cf.h_path,    "%s/H.in", dir);
      sprintf(cf.sigma_path,"%s/Sig.in", dir);
      sprintf(cf.B_path,    "%s/B.out", dir);
      sprintf(cf.V_path,    "%s/V.out", dir);
      /*sprintf(cf.B_path,    "/dev/null");*/
      /*sprintf(cf.V_path,    "/dev/null");*/

      if ( var == 'e' )
          fgls_eigen( 
                  cf.n, cf.p, cf.m, cf.t, cf.wXL, cf.wXR,
                  cf.x_b, cf.y_b, cf.NUM_COMPUTE_THREADS,
                  cf.Phi_path, cf.h_path, cf.sigma_path, 
                  cf.XL_path, cf.XR_path, cf.Y_path, 
                  cf.B_path, cf.V_path
          );
#ifdef FGLS_WITH_GPU
      else if ( var == 'g' )
          fgls_chol_gpu( 
                  cf.n, cf.p, cf.m, cf.t, cf.wXL, cf.wXR,
                  cf.x_b, cf.y_b, cf.NUM_COMPUTE_THREADS,
                  cf.Phi_path, cf.h_path, cf.sigma_path, 
                  cf.XL_path, cf.XR_path, cf.Y_path, 
                  cf.B_path, cf.V_path
          );
#endif
      else
          fgls_chol( 
                  cf.n, cf.p, cf.m, cf.t, cf.wXL, cf.wXR,
                  cf.x_b, cf.y_b, cf.NUM_COMPUTE_THREADS,
                  cf.Phi_path, cf.h_path, cf.sigma_path, 
                  cf.XL_path, cf.XR_path, cf.Y_path, 
                  cf.B_path, cf.V_path
          );
  }

  gettimeofday(&end, NULL);
  long total = get_diff_ms(&start, &end) / nrep; // float?
  printf("%ld\n", total);

#ifdef DEBUG
  printf("Checking B\n");
  /*sleep(3);*/

  char str_buf1[STR_BUFFER_SIZE];
  char str_buf2[STR_BUFFER_SIZE];
  double *mine, *expected;
  FILE *mine_f, *expected_f;
  double max;

  sprintf(str_buf1, "%s/B.out", dir);
  sprintf(str_buf2, "%s/B_exp.out", dir);
  mine_f = fopen(str_buf1, "rb");
  expected_f = fopen(str_buf2, "rb");
  if(mine_f && expected_f) {
    mine = (double*) malloc (cf.m * cf.t * cf.p * sizeof(double));
    expected  = (double*) malloc (cf.m * cf.t * cf.p * sizeof(double));

    sync_read(mine, mine_f, cf.p * cf.m * cf.t, 0);
    sync_read(expected, expected_f,  cf.p * cf.m * cf.t, 0);
    max = compare(mine, expected, cf.m*cf.p*cf.t);
    printf("Max elemental diff: %.16e\n", max);

    free( mine );
    free( expected );
  } else {
    printf("ERROR opening %s and/or %s\n", str_buf1, str_buf2);
  }

  printf("Checking V\n");
  sprintf(str_buf1, "%s/V.out", dir);
  sprintf(str_buf2, "%s/V_exp.out", dir);
  mine_f = fopen(str_buf2, "rb");
  expected_f = fopen(str_buf1, "rb");

  if(mine_f && expected_f) {
    mine = (double *) malloc (cf.m * cf.t * cf.p * cf.p * sizeof(double));
    expected = (double *) malloc (cf.m * cf.t * cf.p * cf.p * sizeof(double));

    sync_read(mine, mine_f, cf.p * cf.p * cf.m * cf.t, 0);
    sync_read(expected,  expected_f,  cf.p * cf.p * cf.m * cf.t, 0);
    max = compare(mine, expected, cf.m * cf.p * cf.p * cf.t);
    printf("Max elemental diff: %.16e\n", max);

    free( mine );
    free( expected );
  } else {
    printf("ERROR opening %s and/or %s\n", str_buf1, str_buf2);
    return -1;
  }

  printf("Done checking\n");
#endif // DEBUG

  return 0;
}

double compare(double* a, double* b, int size) {
  double maxerr = 0.0, relerr = 0.0;
  int i, maxerr_i;
  for (i = 0; i < size; i++) {
    relerr = fabs((a[i]-b[i])/b[i]);
/*    if(relerr > 1e-15)
    {
        printf("Difference at %d: %e [(%f - %f)/%f]\n", i, relerr, a[i], b[i], b[i]);
        return -1;
    }*/
    if(relerr > maxerr) {
      maxerr = relerr;
      maxerr_i = i;
    }
  }
  if(maxerr > 1e-15)
    printf("max elemental relative error of %g at %d. Maximum absolute error there: %g (%.15e-%.15e)\n", maxerr, maxerr_i, a[maxerr_i]-b[maxerr_i], a[maxerr_i], b[maxerr_i]);
  return maxerr;
}

