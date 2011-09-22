#include "fgls.h"
#include "io.h"
#include "fgls_eigen.h"

#include <malloc.h>
#include <stdio.h>

int main( int argc, char *argv[] ) 
{
  char *dir;
  FGLS_eigen_t *cf = &FGLS_eigen_config;

  if (argc != 2) {
    printf("usage: %s <data_set_directory>\n\tnote: data_set_directory must contain files X.in, Y.in, Phi.in, and H.in (B_exp.out if checking answer)\n", argv[0]);
    return -1;
  }

  dir = argv[1];
  printf("Please enter parameters\n");
  printf("\tm: ");
  scanf("%d", &cf->m);
  printf("\tt: ");
  scanf("%d", &cf->t);
  printf("\tm blocksize: ");
  scanf("%d", &cf->x_b);
  // Still needed for OOC gemm
  /*printf("\tt blocksize: ");*/
  /*scanf("%d", &cf->y_b);*/
  printf("\tn: ");
  scanf("%d", &cf->n);
  printf("\tp: ");
  scanf("%d", &cf->p);
  printf("\tcompute threads: ");
  scanf("%d", &cf->NUM_COMPUTE_THREADS);

  /*in.m_indexed = (int) ((double) in.m/in.x_b+.5);*/
  /*in.t_indexed = (int) ((double) in.t/in.y_b+.5);*/
  /*in.NUM_BUFFERS_PER_THREAD = 2;*/

#if TIMING
  cf->time = (timing*)malloc(sizeof(timing));
  struct timeval start, end;
  gettimeofday(&start, NULL);
#endif // TIMING

  sprintf(cf->X_path,   "%s/X.in", dir);
  sprintf(cf->ZtX_path, "%s/X.tmp", dir);
  sprintf(cf->Y_path,   "%s/Y.in", dir);
  sprintf(cf->ZtY_path, "%s/Y.tmp", dir);
  sprintf(cf->Phi_path, "%s/Phi.in", dir);
  sprintf(cf->h_path,   "%s/H.in", dir);
  sprintf(cf->B_path,   "%s/B.out", dir);

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
  read(b_mine, b_mine_f, cf->p*cf->m*cf->t, 0);
  read(b_exp, b_exp_f, cf->p*cf->m*cf->t, 0);
  printf("out[0]: %12e\n", b_mine[0]);
  printf("exp[0]: %12e\n", b_exp[0]);
  printf("out[1]: %12e\n", b_mine[1]);
  printf("exp[1]: %12e\n", b_exp[1]);
  /*printf("out[2]: %12e\n", b_mine[2]);*/
  /*printf("exp[2]: %12e\n", b_exp[2]);*/
  /*printf("out[3]: %12e\n", b_mine[3]);*/
  /*printf("exp[3]: %12e\n", b_exp[3]);*/
  /*printf("out[4]: %12e\n", b_mine[4]);*/
  /*printf("exp[4]: %12e\n", b_exp[4]);*/
  /*printf("out[400]: %12e\n", b_mine[400]);*/
  /*printf("exp[400]: %12e\n", b_exp[400]);*/
  printf("out[3000]: %12e\n", b_mine[3000]);
  printf("exp[3000]: %12e\n", b_exp[3000]);
  printf("out[6000]: %12e\n", b_mine[6000]);
  printf("exp[6000]: %12e\n", b_exp[6000]);
  printf("out[9000]: %12e\n", b_mine[9000]);
  printf("exp[9000]: %12e\n", b_exp[9000]);
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
