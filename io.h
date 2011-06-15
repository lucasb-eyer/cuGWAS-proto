#include <stdio.h>


int write_float(float *buffer, FILE* file, int stride, int num_indices, int start_index);
int write_double(double *buffer, FILE* file, int stride, int num_indices, int start_index);

int read_float(float *buffer, FILE* file, int stride, int num_indices, int start_index);
int read_double(double *buffer, FILE* file, int stride, int num_indices, int start_index);

void print_buffer(double *buffer, int items);
