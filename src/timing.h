#ifndef TIMING_H
#define TIMING_H

#include <stdio.h>

#include <sys/time.h>
#include <time.h>

long get_diff_ms(struct timeval *s, struct timeval *e);
long get_diff_us(struct timeval *s, struct timeval *e);
long get_diff_ns(struct timespec *s, struct timespec *e);

#endif
