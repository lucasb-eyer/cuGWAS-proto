#include "fgls.h"
#include <stdio.h>
#include <stddef.h>

void write(double *buffer, FILE* file, size_t size, size_t location);

void read(double *buffer, FILE* file, size_t size, size_t start);

void print_buffer(double *buffer, int items);
