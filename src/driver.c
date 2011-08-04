#include "fgls.h"
#include "io.h"
#include "fgls_eigen.h"

#include <malloc.h>
#include <stdio.h>

int main(int argc, char* argv[]) {
  char trav[1];
  problem_args in;
  char *dir;
  int eigen = 0;

  if (argc != 2) {
    printf("usage: %s <data_set_directory>\n\tnote: data_set_directory must contain files X.in, Y.in, Phi.in, and H.in (B_exp.out if checking answer)\n", argv[0]);
    return -1;
  }

  dir = argv[1];
  printf("Please enter parameters\n");
  printf("\tm: ");
  scanf("%d", &in.m);
  printf("\tt: ");
  scanf("%d", &in.t);
  printf("\tm blocksize: ");
  scanf("%d", &in.x_b);
  // Still needed for OOC gemm
  printf("\tt blocksize: ");
  scanf("%d", &in.y_b);
  printf("\tn: ");
  scanf("%d", &in.n);
  printf("\tp: ");
  scanf("%d", &in.p);

  in.m_indexed = (int) ((double) in.m/in.x_b+.5);
  in.t_indexed = (int) ((double) in.t/in.y_b+.5);

#if TIMING
  in.time = (timing*)malloc(sizeof(timing));
  struct timeval start, end;
  gettimeofday(&start, NULL);
#endif // TIMING

  fgls_eigen(dir, &in);

#if TIMING
  gettimeofday(&end, NULL);
  long total = get_diff_ms(&start, &end);;
#endif // TIMING
#if DEBUG
  char str_buf[STR_BUFFER_SIZE];
  double *b_mine, *b_exp;
  FILE *b_mine_f, *b_exp_f;
  b_mine = (double*)malloc(in.m*in.t*in.p*sizeof(double));
  b_exp = (double*)malloc(in.m*in.t*in.p*sizeof(double));

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
  read(b_mine, b_mine_f, in.p*in.m*in.t, 0);
  read(b_exp, b_exp_f, in.p*in.m*in.t, 0);
  max = compare(b_mine, b_exp, in.m*in.t*in.p);
  printf("Max elemental diff: %lf\n", max);
#endif // DEBUG

#if TIMING
  if (!in.time) {
    return 0;
  }
  printf("%c timing_info [ m t ] = [ total io comp io_wait comp_wait ];\n", '%');
  printf("timing_info [ %d %d ] = [ %ld %ld %ld %ld %ld ];\n", in.x_b, in.y_b, total, in.time->io_time, in.time->compute_time, 
         in.time->io_mutex_wait_time, in.time->comp_mutex_wait_time);
#endif // TIMING
  return 0;
}
