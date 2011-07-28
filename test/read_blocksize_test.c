#include "read_blocksize_test.h"

#include "src/io.h"

#include <stdio.h>
#include <sys/time.h>
#include <time.h>

int main(int argc, char* argv[]) {
  test_params i = parse_args(argc, argv);
  read_blocksize_test(&i);
  cleanup_test_params(&i);
}

void read_blocksize_test(test_params* in) {
  //////// declarations ////////
  fprintf(stdout, "%c start read_blocksize_test\n", '%');
  char *colors = "brkgmcbrkg";
  char *ticks  = "o+*xso+*xs";

  struct timeval start_t, end_t;
  double* buf;
  int i, p;
  long mtime, mtime_max = 0;

  //////// run the test //////// 
  for (p = 1, i = in->start; i < in->end; i += in->inc, p++) {
    buf = (double*) malloc( in->size * sizeof(double));
    fprintf(in->output_file, "read_blocksize_test_( %d, 1:2 ) = [ %d  ", p, i/1000);

    int num_blocks = in->size/i;
    int last_block_size = in->size%i;

    gettimeofday(&start_t, NULL);
    read(buf, in->input_file, num_blocks*i, 0);

    // handle corner case of blocksize/file size mismatch
    if(last_block_size != 0) 
      read(buf, in->input_file, last_block_size, num_blocks*i);    

    gettimeofday(&end_t, NULL);
    mtime = get_diff_ms(&start_t, &end_t);

    fprintf(in->output_file, "%lu ];\n", mtime);
    if (mtime > mtime_max) 
      mtime_max = mtime;
    free(buf);
  }

  //////// print matlab formatting info ////////
  fprintf(in->output_file, "figure;\n");
  fprintf(in->output_file, "hold on;\n");

  fprintf(in->output_file, "plot( read_blocksize_test_( :,1 ), read_blocksize_test_( :, 2 ), '%c:%c' ); \n",
	  colors[ 0 ], ticks[ 0 ]);

  fprintf(in->output_file, "xlabel( 'blocksize i' );\n");
  fprintf(in->output_file, "ylabel( 'milliseconds' );\n");
  fprintf(in->output_file, "axis( [ %ld %ld 0 %lu ] ); \n", in->start, in->end, mtime_max);
  fprintf(in->output_file, "title( 'FLGS Read Blocksize Test' );");
  fprintf(in->output_file, "print -depsc flgs_read_blocksize_test.eps;\n");
  fprintf(in->output_file, "hold off;\n");
}
