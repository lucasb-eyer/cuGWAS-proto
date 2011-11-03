#ifndef IO_H
#define IO_H

#include <stdio.h>

void sync_write(double *buffer, FILE* file, size_t size, size_t start);

void sync_read(double *buffer, FILE* file, size_t size, size_t start);

#endif
