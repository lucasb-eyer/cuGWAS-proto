#include "test_framework.h"
#include "read_tests.h"
#include "write_tests.h"

#include <stdio.h>

void run_test(int test) {
  switch (test) {
  case READ_TEST:
    read_test();
    break;
  case READ_BLOCKSIZE_TEST:
    read_blocksize_test();
    break;
  case WRITE_TEST:
    write_test();
    break;
  case WRITE_BLOCKSIZE_TEST:
    write_blocksize_test();
    break;
  default:
    break;
  }
}

void run_all_tests() {
  run_test(READ_TEST);
  run_test(READ_BLOCKSIZE_TEST);
  run_test(WRITE_TEST);
  run_test(WRITE_BLOCKSIZE_TEST);
}
