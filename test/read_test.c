#include "src/io.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void read_test() {
  char *colors = "brkgmcbrkg";
  char *ticks  = "o+*xso+*xs";

  int i, p, n_repeats, start, end;

  FILE* fp;

  double max_gflops=6.0;

  time_t dtime_first, dtime_last, dtime_max;

  double
    gflops,
    diff;

  fp = fopen("test/input", "rb");
  if (!fp) {
    fprintf(stderr, "file open error in read_test");
    exit(-1)
  }

  fprintf( stdout, "%c number of repeats: ", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf(stdout, "%c enter start and end problem sizes ", '%');
  scanf("%d%d", &start, &end);
  fprintf(stdout, "%c %d %d\n", '%', start, end);
  fprintf(stdout, "\n");

  for (p = 0, i=start; i < end; i*=2, p++) {
    fprintf(stdout, "read_test( %d, 1:3 ) = [ %d  ", p, i);
    fflush(stdout);
    ctime(dtime_first);
    read_double(buf, fp, 1, i, 0);
    ctime(dtime_last);
    fprintf(stdout, "%6.3lf", dtime_last - dtime_first);
    fprintf(stdout, " ]; \n");
    fflush(stdout);
    if (dtime_last - dtime_first > dtime_max) 
      dtime_max = dtime_last - dtime_first;
  }

  fprintf(stdout, "figure;\n");
  fprintf(stdout, "hold on;\n");

  fprintf(stdout, "plot( read_test( :,1 ), read_test( :, 2 ), '%c:%c' ); \n",
	  colors[ 0 ], ticks[ 0 ]);

  fprintf(stdout, "xlabel( 'input size i' );\n");
  fprintf(stdout, "ylabel( 'GFLOPS/sec.' );\n");
  fprintf(stdout, "axis( [ %d %d 0 %d ] ); \n", start, end, dtime_max);
  fprintf(stdout, "title( 'FLGS Read Test' );");
  fprintf(stdout, "print -depsc flgs_read_test.eps;\n");
  fprintf(stdout, "hold off;\n");
  fflush(stdout);
}
