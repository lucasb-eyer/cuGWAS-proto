#include "test_framework.h"

#include <stdio.h>
#include <stdlib.h>

#define BUF_SIZE 256

void randomize(double* buf, int size) {
  int i;
  for( i = 0; i < size; i++) {
    buf[i] = 1000*rand();
  }
}



// <repeats> <size> <start> <end> <inc>
test_params parse_args(int argc, char* argv[]) {
  test_params out;

  //////// arg parsing ////////
  if (argc != 6) {
      printf("usage: %s <repeats> <size> <start> <end> <inc>\n", argv[0]);
      exit(-1);
  }

  out.repeats = atoi(argv[1]);
  out.size = atoll(argv[2]);
  out.start = atoll(argv[3]);
  out.end = atoll(argv[4]);
  out.inc = atoll(argv[5]);

  //////// get io files ////////
  char buffer[BUF_SIZE];
  printf("%c enter input file: ", '%');
  fgets(buffer, BUF_SIZE, stdin);
  out.input_file = fopen(buffer, "wb+");
  if (!out.input_file) {
    fprintf(stderr, "error opening %s\n", buffer);
    exit(-1);
  }

  printf("%c enter output file: ", '%');
  fgets(buffer, BUF_SIZE, stdin);
  out.output_file = fopen(buffer, "w");
  if (!out.output_file) {
    fprintf(stderr, "error opening %s\n", buffer);
    exit(-1);
  }

  return out;
}

void cleanup_test_params(test_params* i) {
  if (i) {
    if (i->input_file) {
      fclose(i->input_file);
      i->input_file = 0;
    }
    if (i->output_file) {
      fclose(i->output_file);
      i->output_file = 0;
    }
  }
}
