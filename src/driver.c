#include "fgls.h"
#include "io.h"
#include "t_traversal_eigen.h"
#include "t_traversal_chol.h"
#include "m_traversal_eigen.h"
#include "m_traversal_chol.h"
#include "test/test_framework.h"
#include <malloc.h>
#include <stdio.h>

#define M_MAX 1000
#define T_MAX 15

int main(int argc, char* argv[]) {
  char trav[1];
  problem_args in;
  char *x, *y, *phi, *b;
  int eigen = 0;

  if (argc != 7) {
    printf("usage: %s <eigen|chol> <m|t> <x-in-file> <y-in-file> <phi-in-file> <b-out-file>\n", argv[0]);
    return -1;
  }

  eigen = argv[1][0] == 'e';
  trav[0] = argv[2][0];
  x = argv[3];
  y = argv[4];
  phi = argv[5];
  b = argv[6];
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
  printf("\th: ");
  scanf("%lf", &in.h);


  in.m_indexed = (int) ((double) in.m/in.x_b);
  in.t_indexed = (int) ((double) in.t/in.y_b);

#if TIMING
  in.time = (timing*)malloc(sizeof(timing));
#endif // TIMING

  if (trav[0] == 'm') {
    if (eigen) {
      m_traversal_eigen(x, y, phi, b, &in);
    } else {
      m_traversal_chol(x, y, phi, b, &in);
    }
  } else {
    if (eigen) {
      t_traversal_eigen(x, y, phi, b, &in);
    } else {
      t_traversal_chol(x, y, phi, b, &in);
    }
  }
  printf("here\n");
#if TIMING
  double *b_mine, *b_exp;
  FILE *b_mine_f, *b_exp_f;
  b_mine = (double*)malloc(in.p*sizeof(double));
  b_exp = (double*)malloc(in.p*sizeof(double));
  b_mine_f = fopen("b.out", "rb");
  if(!b_mine_f) {
    printf("ERROR opening my b\n");
    return -1;
  }
  b_exp_f = fopen("B_exp.out", "rb");
  if(!b_exp_f) {
    printf("ERROR opening exp b\n");
    return -1;
  }
  double temp, max = 0.0;
  int j, i, my_index, exp_index;
  for (i = 0; i < in.m; i++) {
    for (j = 0; j < in.t; j++) {
      my_index = in.m*j + i;
      exp_index = M_MAX*j + i;
      read(b_mine, b_mine_f, in.p, my_index);
      read(b_exp, b_exp_f, in.p, exp_index);
      printf("my: %d exp: %d\n", my_index, exp_index);
      temp = compare(b_mine, b_exp, in.p);
      if( temp > max )
        max = temp;
    }  
  }
  printf("Max elemental diff: %lf\n", max);
#endif // DEBUG

#if TIMING
  if (!in.time) {
    return 0;
  }
  printf("timing results:\n");
  printf("\ttime in io: %ld\n", in.time->io_time);
  printf("\ttime in compute: %ld\n", in.time->compute_time);
  printf("\ttime in io_mutex: %ld\n", in.time->io_mutex_wait_time);
  printf("\ttime in compute_mutex: %ld\n", in.time->comp_mutex_wait_time);

#endif // TIMING
  return 0;
}


