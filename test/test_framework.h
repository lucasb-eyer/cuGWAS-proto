#ifndef TEST_FRAMEWORK_H
#define TEST_FRAMEWORK_H

#define READ_TEST              0
#define READ_BLOCKSIZE_TEST    1         
#define WRITE_TEST             2
#define WRITE_BLOCKSIZE_TEST  3
#define X_TRAVERSAL_TEST       4
#define Y_TRAVERSAL_TEST       5

typedef struct test_timings_t {
  unsigned int comp_mutex_wait;
  unsigned int compute_time;
  unsigned int io_mutex_wait;
  unsigned int io_time;
} test_timings;

void run_test(int test);
void run_all_tests();


#endif // TEST_FRAMEWORK_H 
