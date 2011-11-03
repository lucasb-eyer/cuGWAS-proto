#ifndef TIMING_H
#define TIMING_H

#include <stdio.h>

#include <sys/time.h>
#include <time.h>

#if TIMING
#define BEGIN_TIMING() \
	gettimeofday(&start, NULL);
#define END_TIMING(var)     \
	gettimeofday(&end, NULL);   \
    var += get_diff_ms(&start, &end);
#define DEF_TIMING() \
	struct timeval start, end;
#else // TIMING
#define BEGIN_TIMING()
#define END_TIMING(var)
#define DEF_TIMING()
#endif // TIMING

typedef struct timing_t {
  long compute_time;
  long io_time;
  long comp_mutex_wait_time;
  long io_mutex_wait_time;
} timing;

long get_diff_ms(struct timeval *s, struct timeval *e);
long get_diff_us(struct timeval *s, struct timeval *e);
long get_diff_ns(struct timespec *s, struct timespec *e);

#endif
