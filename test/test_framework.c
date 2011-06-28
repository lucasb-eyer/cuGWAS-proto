#include "test_framework.h"
#include "read_tests.h"
#include "write_tests.h"

#include <stdio.h>

void run_test(int test, int repeats, int size, int start, int end, int inc) {
  switch (test) {
  case READ_TEST:
    read_test(repeats, start, end, inc);
    break;
  case READ_BLOCKSIZE_TEST:
    read_blocksize_test(repeats, size, start, end, inc);
    break;
  case WRITE_TEST:
    write_test(repeats, start, end, inc);
    break;
  case WRITE_BLOCKSIZE_TEST:
    write_blocksize_test(repeats, size, start, end, inc);
    break;
  default:
    break;
  }
}

void run_all_tests(int repeats, int size, int start, int end, int inc) {
  run_test(READ_TEST, repeats, size, start, end, inc);
  run_test(READ_BLOCKSIZE_TEST, repeats, size, start, end, inc);
  run_test(WRITE_TEST, repeats, size, start, end, inc);
  run_test(WRITE_BLOCKSIZE_TEST, repeats, size, start, end, inc);
}
