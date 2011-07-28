#include "write_test.h"

#include "src/io.h"

#include <stdio.h>
#include <sys/time.h>
#include <time.h>

int main(int argc, char* argv[]) {
  test_params i = parse_args(argc, argv);
  write_test(&i);
  cleanup_test_params(&i);
}

void write_test(test_params* in) {
  //////// declarations ////////
  fprintf(stdout, "%c start write_test\n", '%');
  char *colors = "brkgmcbrkg";
  char *ticks  = "o+*xso+*xs";

  struct timeval start_t, end_t;
  double* buf;
  int i, p;
  long mtime, mtime_max = 0;

  //////// run the test //////// 
  buf = (double*) malloc(in->end * sizeof(double));
  randomize(buf, in->end );
  for (p = 1, i = in->start; i < in->end; i += in->inc, p++) {
    fprintf(in->output_file, "write_test_( %d, 1:2 ) = [ %d  ", p, i);

    gettimeofday(&start_t, NULL);
    write(buf, in->input_file, in->size, i);
    gettimeofday(&end_t, NULL);
    mtime = get_diff_ms(&start_t, &end_t);  

    fprintf(in->output_file, "%lu ]; \n", mtime);
    if (mtime > mtime_max) 
      mtime_max = mtime;
  }
  free(buf);

  //////// print matlab formatting info ////////
  fprintf(in->output_file, "figure;\n");
  fprintf(in->output_file, "hold on;\n");

  fprintf(in->output_file, "plot( write_test_( :,1 ), write_test_( :, 2 ), '%c:%c' ); \n",
	  colors[ 0 ], ticks[ 0 ]);

  fprintf(in->output_file, "xlabel( 'input size i' );\n");
  fprintf(in->output_file, "ylabel( 'milliseconds' );\n");
  fprintf(in->output_file, "axis( [ %ld %ld 0 %lu ] ); \n", in->start, in->end, mtime_max);
  fprintf(in->output_file, "title( 'FLGS Write Test' );");
  fprintf(in->output_file, "print -depsc flgs_write_test.eps;\n");
  fprintf(in->output_file, "hold off;\n");
}
