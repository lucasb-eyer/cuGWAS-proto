#include <stdio.h>
#include "fgls.h"

int write(double *buffer, FILE* file, int size, int location);

int read(double *buffer, FILE* file, int stride, int num_indices, int start_index);

void print_buffer(double *buffer, int items);
