#include "fgls.h"
#include <stdio.h>
#include <stddef.h>

void my_write(double *buffer, FILE* file, size_t size, size_t location);

void my_read(double *buffer, FILE* file, size_t size, size_t start);

/*void print_buffer(double *buffer, int items);*/
