#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sys/time.h>
#include <time.h>

#include "io.h"
#include "timing.h"

long get_diff_ms(struct timeval *start_time, struct timeval *end_time) {
  long seconds = end_time->tv_sec - start_time->tv_sec;
  long useconds = end_time->tv_usec - start_time->tv_usec;
  return ((seconds) * 1000 + useconds/1000.0) +0.5;
}

long get_diff_us(struct timeval *start_time, struct timeval *end_time) {
  long seconds = end_time->tv_sec - start_time->tv_sec;
  long useconds = end_time->tv_usec - start_time->tv_usec;
  return ((seconds) * 1000000 + useconds);
}
long get_diff_ns(struct timespec *start, struct timespec *end)
{
	long seconds, nseconds;
	long nsecs_in_sec = 1e9;
	if ( (end->tv_nsec - start->tv_nsec) < 0 ) 
	{
        seconds = end->tv_sec - start->tv_sec - 1;
		nseconds = nsecs_in_sec - (end->tv_nsec - start->tv_nsec);
		return nsecs_in_sec * seconds + nseconds;
	} 
	else 
	{
        seconds = end->tv_sec - start->tv_sec;
		nseconds = end->tv_nsec - start->tv_nsec;
		return nsecs_in_sec * seconds + nseconds;
	}
}
