#include "fgls.h"
#include "test/test_framework.h"
#include "test/io_startup_test.h"
#include <stdio.h>

int main(int argc, char* argv[]) {
  /*  int n_repeats, size, start, end, inc;

  fprintf( stdout, "%c number of repeats: ", '%' );
  scanf( "%d", &n_repeats );

  fprintf(stdout, "%c enter array size, start, end, and inc blocksize(in KB): ", '%');
  scanf("%d%d%d%d", &size, &start, &end, &inc);

  // adjust for kilobytes
  start*=1000;
  end*=1000;
  size*=1000;
  inc*=1000;

//  run_all_tests(n_repeats, size, start, end, inc);*/
  io_startup_test();

  return 0;
}


