#include "fgls.h"
#include "io.h"
#include "t_traversal_eigen.h"
#include "t_traversal_chol.h"
#include "m_traversal_eigen.h"
#include "m_traversal_chol.h"
//#include "test/test_framework.h"
#include <malloc.h>
#include <stdio.h>
/*
int main(int argc, char* argv[]) {
  char trav[1];
  problem_args in;
  char *x, *y, *phi, *b, *h;
  int eigen = 0;

  if (argc != 8) {
    printf("usage: %s <eigen|chol> <m|t> <x-in-file> <y-in-file> <phi-in-file> <h-in-file> <b-out-file>\n", argv[0]);
    return -1;
  }

  eigen = argv[1][0] == 'e';
  trav[0] = argv[2][0];
  x = argv[3];
  y = argv[4];
  phi = argv[5];
  h = argv[6];
  b = argv[7];
  printf("Please enter parameters\n");
  printf("\tm: ");
  scanf("%d", &in.m);
  printf("\tt: ");
  scanf("%d", &in.t);
  printf("\tm blocksize: ");
  scanf("%d", &in.x_b);
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


  if (trav[0] == 'm') {
    if (eigen) {
      m_traversal_eigen(x, y, phi, h, b, &in);
    } else {
      m_traversal_chol(x, y, phi, h, b, &in);
    }
  } else {
    if (eigen) {
      t_traversal_eigen(x, y, phi, h, b, &in);
    } else {
      t_traversal_chol(x, y, phi, h, b, &in);
    }
  }
#if TIMING
  gettimeofday(&end, NULL);
  long total = get_diff_ms(&start, &end);;
#endif // TIMING
#if DEBUG
  double *b_mine, *b_exp;
  FILE *b_mine_f, *b_exp_f;
  b_mine = (double*)malloc(in.m*in.t*in.p*sizeof(double));
  b_exp = (double*)malloc(in.m*in.t*in.p*sizeof(double));
  b_mine_f = fopen("/home/rt203005/rt203005_FGLS/trunk/b.out", "rb");
  if(!b_mine_f) {
    printf("ERROR opening my b\n");
    return -1;
  }
  b_exp_f = fopen("/home/rt203005/rt203005_FGLS/trunk/B_exp.out", "rb");
  if(!b_exp_f) {
    printf("ERROR opening exp b\n");
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


*/
