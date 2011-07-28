#ifndef TEST_FRAMEWORK_H
#define TEST_FRAMEWORK_H

#include <stddef.h>
#include <stdio.h>

typedef struct test_timings_t {
  unsigned int comp_mutex_wait;
  unsigned int compute_time;
  unsigned int io_mutex_wait;
  unsigned int io_time;
} test_timings;

typedef struct test_params_t {
  int repeats;
  size_t size;
  size_t start;
  size_t end;
  size_t inc;
  FILE* input_file;
  FILE* output_file;
} test_params;

//////// matrix operations ////////
void randomize(double* buf, int size);

//////// parameter functions ////////
test_params parse_args(int argc, char* argv[]);
void cleanup_test_params(test_params* i);

#endif // TEST_FRAMEWORK_H 
