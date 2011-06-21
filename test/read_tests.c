#include "read_tests.h"

#include "src/io.h"

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void read_test() {
  fprintf( stdout, "%c start read_test\n", '%' );
  char *colors = "brkgmcbrkg";
  char *ticks  = "o+*xso+*xs";

  double* buf;
 
  int i, p, n_repeats, start, end, inc;

  FILE* fp;
  FILE* out;
  double max_gflops=6.0;

  struct timeval start_time, end_time;
  long mtime, seconds, useconds, mtime_max = 0;

  double
    gflops,
    diff;

  fp = fopen("test/input", "rb");
  if (!fp) {
    fprintf(stderr, "file open error in read_test");
    exit(-1);
  }
  out = fopen("read_test.m", "w");
  if (!out) {
    fprintf(stderr, "file open error in read_test");
    exit(-1);
  }

  fprintf( stdout, "%c number of repeats: ", '%' );
  scanf( "%d", &n_repeats );
  fprintf( out, "%c %d\n", '%', n_repeats );

  fprintf(stdout, "%c enter start, end, and inc (in KB): ", '%');
  scanf("%d%d%d", &start, &end, &inc);
  fprintf(out, "%c %d %d %d\n", '%', start, end, inc);
  fprintf(out, "\n");

  // adjust for kilobytes
  start*=1000;
  end*=1000;
  inc*=1000;

  for (p = 1, i=start; i < end; i+=inc, p++) {
    buf = (double*) malloc( i * sizeof(double));
    fprintf(out, "read_test( %d, 1:2 ) = [ %d  ", p, i);
    fflush(out);
    gettimeofday(&start_time, NULL);
    read_double(buf, fp, 1, i, 0);
    gettimeofday(&end_time, NULL);
    seconds = end_time.tv_sec - start_time.tv_sec;
    useconds = end_time.tv_usec - start_time.tv_usec;
    mtime = ((seconds) * 1000 + useconds/1000.0) +0.5;
    fprintf(out, "%lu", mtime);
    fprintf(out, " ]; \n");
    fflush(out);
    if (mtime > mtime_max) 
      mtime_max = mtime;
    free(buf);
  }

  fprintf(out, "figure;\n");
  fprintf(out, "hold on;\n");

  fprintf(out, "plot( read_test( :,1 ), read_test( :, 2 ), '%c:%c' ); \n",
	  colors[ 0 ], ticks[ 0 ]);

  fprintf(out, "xlabel( 'input size i' );\n");
  fprintf(out, "ylabel( 'GFLOPS/sec.' );\n");
  fprintf(out, "axis( [ %d %d 0 %lu ] ); \n", start, end, mtime_max);
  fprintf(out, "title( 'FLGS Read Test' );");
  fprintf(out, "print -depsc flgs_read_test.eps;\n");
  fprintf(out, "hold off;\n");
  fflush(out);

  fclose(fp);
  fclose(out);
}


void read_blocksize_test() {
  fprintf( stdout, "%c start read_blocksize_test\n", '%' );

  char *colors = "brkgmcbrkg";
  char *ticks  = "o+*xso+*xs";

  double* buf;
 
  int i, p, n_repeats, start, end, size, inc;

  FILE* fp;
  FILE* out;
  double max_gflops=6.0;

  struct timeval start_time, end_time;
  long mtime, seconds, useconds, mtime_max = 0;

  double
    gflops,
    diff;

  fp = fopen("test/input", "rb");
  if (!fp) {
    fprintf(stderr, "file open error in read_test");
    exit(-1);
  }

  out = fopen("read_blocksize_test", "w");
  if (!out) {
    fprintf(stderr, "file open error in read_test");
    exit(-1);
  }

  fprintf( stdout, "%c number of repeats: ", '%' );
  scanf( "%d", &n_repeats );
  fprintf( out, "%c %d\n", '%', n_repeats );

  fprintf(stdout, "%c enter array size, start, end, and inc blocksize(in KB): ", '%');
  scanf("%d%d%d%d", &size, &start, &end, &inc);
  fprintf(out, "%c %d %d %d %d\n", '%', size, start, end, inc);
  fprintf(out, "\n");

  // adjust for kilobytes
  start*=1000;
  end*=1000;
  size*=1000;
  inc*=1000;

  for (p = 1, i=start; i < end; i+=inc, p++) {
    buf = (double*) malloc( size * sizeof(double));
    fprintf(out, "read_blocksize_test( %d, 1:2 ) = [ %d  ", p, i/1000);
    fflush(out);
    int num_blocks = size/i;
    int last_block_size = size%i;
    gettimeofday(&start_time, NULL);
    read_double(buf, fp, num_blocks, i, 0);
    if( last_block_size != 0 ) 
      read_double(buf, fp, 1, last_block_size, num_blocks*i);
    gettimeofday(&end_time, NULL);
    seconds = end_time.tv_sec - start_time.tv_sec;
    useconds = end_time.tv_usec - start_time.tv_usec;
    mtime = ((seconds) * 1000 + useconds/1000.0) +0.5;
    fprintf(out, "%lu", mtime);
    fprintf(out, " ]; \n");
    fflush(out);
    if (mtime > mtime_max) 
      mtime_max = mtime;
    free(buf);
  }

  fprintf(out, "figure;\n");
  fprintf(out, "hold on;\n");

  fprintf(out, "plot( read_blocksize_test( :,1 ), read_blocksize_test( :, 2 ), '%c:%c' ); \n",
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
