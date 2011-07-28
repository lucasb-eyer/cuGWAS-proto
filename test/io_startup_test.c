#include "test_framework.h"

#include <stdlib.h>
#include <stdio.h>

void write_io(char* f, size_t size) {
  FILE* fp;
  int i;
  double *buf;
  fp = fopen(f, "wb");
  if (!fp) 
    printf("OH NOES. IO FAIL WITH %s\n", f);
  buf = (double*) malloc(size * sizeof(double));
  for (i = 0; i < size; i++) {
    buf[i] = rand() *100;
  }
  fwrite(buf, size*sizeof(double), 1, fp);
  fclose(fp);
  free(buf);
}
void io_startup_test() {

  write_io("test/input", 1000000000);

  // read 1kb block at every 1mb in a file that is 1gb
  run_test(READ_TEST, 1, 1000, 0, 1000000000, 10000000); 

  // write 1kb block at every 1mb in a file that is 1gb
  //  run_test(WRITE_BLOCKSIZE_TEST, 1, 1000, 0, 1000000000, 1000000); 
}


