#include <stdio.h>
#include "fgls.h"

int write(double *buffer, FILE* file, int size, int location);

int read(double *buffer, FILE* file, int size, int start);

void print_buffer(double *buffer, int items);
