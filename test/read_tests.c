#include "read_tests.h"

#include "src/io.h"
#include "src/fgls.h"

#include <stdio.h>
#include <stdlib.h>


void read_blocksize_test(int n_repeats, size_t size, size_t start, size_t end, size_t inc) {
  fprintf( stdout, "%c start read_blocksize_test\n", '%' );

  char *colors = "brkgmcbrkg";
  char *ticks  = "o+*xso+*xs";

  double* buf;
 
  int i, p;

  FILE* fp;
  FILE* out;
  struct timeval start_time, end_time;
  long mtime, mtime_max = 0;

  fp = fopen("/home/rt203005/rt203005_FGLS/trunk/test/input", "rb");
  if (!fp) {
    fprintf(stderr, "file open error in read_test");
    exit(-1);
  }

  out = fopen("/home/rt203005/rt203005_FGLS/trunk/read_blocksize_test.m", "w");
  if (!out) {
    fprintf(stderr, "file open error in read_test");
    exit(-1);
  }

  for (p = 1, i=start; i < end; i+=inc, p++) {
    buf = (double*) malloc( size * sizeof(double));
    fprintf(out, "read_blocksize_test_( %d, 1:2 ) = [ %d  ", p, i/1000);
    fflush(out);
    int num_blocks = size/i;
    int last_block_size = size%i;
    gettimeofday(&start_time, NULL);
    read(buf, fp, num_blocks*i, 0);
    if( last_block_size != 0 ) 
      read(buf, fp, last_block_size, num_blocks*i);    
    gettimeofday(&end_time, NULL);
    mtime = get_diff_ms(&start_time, &end_time);
    fprintf(out, "%lu", mtime);
    fprintf(out, " ]; \n");
    fflush(out);
    if (mtime > mtime_max) 
      mtime_max = mtime;
    free(buf);
  }

  fprintf(out, "figure;\n");
  fprintf(out, "hold on;\n");

  fprintf(out, "plot( read_blocksize_test_( :,1 ), read_blocksize_test_( :, 2 ), '%c:%c' ); \n",
	  colors[ 0 ], ticks[ 0 ]);

  fprintf(out, "xlabel( 'blocksize i(in kb for size %d)' );\n", size/1000);
  fprintf(out, "ylabel( 'millisec.' );\n");
  fprintf(out, "axis( [ %d %d 0 %lu ] ); \n", start/1000, end/1000, mtime_max);
  fprintf(out, "title( 'FLGS Read Blocksize Test' );");
  fprintf(out, "print -depsc flgs_read_test.eps;\n");
  fprintf(out, "hold off;\n");
  fflush(out);

  fclose(fp);
  fclose(out);
}
