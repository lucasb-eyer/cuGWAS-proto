#include "test_framework.h"

void run_test(int test) {
  switch (test) {
  case READ_TEST:
    read_test();
    break;
  default:
    break;
  }
}
