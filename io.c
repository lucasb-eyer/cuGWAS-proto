#include "io.h"
#include <stdio.h>

int write_general(void* buffer, char* filename, int size, int start) {
  FILE* fp;

  fp = fopen(filename, "wb");
  if(!fp) {
    printf("error: fp null\n");
  }
  fseek(fp, start*sizeof(double), SEEK_SET);
  //just write bytes
  size_t out = fwrite(buffer, sizeof(double), size, fp);
  if(out != size) {
    printf("error: actual write size != proposed write size");
  }
  return fclose(fp);
}
int write_float(float *buffer, char* filename, int length, int num_indices, int start_index) {
  return write_general((void*)buffer, filename, num_indices*length*sizeof(float), start_index*length*sizeof(float));
}
int write_double(double *buffer, char* filename, int length, int num_indices, int start_index) {
  return write_general((void*)buffer, filename, num_indices*length, start_index*length);
}

int read_general(void* buffer, char* filename, int size, int start) {
  FILE* fp;

  fp = fopen(filename, "rb");
  fseek(fp, start*sizeof(double), SEEK_SET);
  //just write bytes
  size_t out = fread(buffer, sizeof(double), size, fp);
  printf("read %x bytes into %x\n", size, buffer);
  if(out != size) {
    printf("error: actual read size != proposed read size");
  }
  return fclose(fp);
}
int read_float(float *buffer, char* filename, int length, int num_indices, int start_index) {
  return read_general((void*)buffer, filename, num_indices*length*sizeof(float), start_index*length*sizeof(float));
}
int read_double(double *buffer, char* filename, int length, int num_indices, int start_index) {
  return read_general((void*)buffer, filename, num_indices*length, start_index*length);
}

void print_buffer(double *buffer, int items) {
  int i;
  for(i = 0; i < items; i++) {
    printf("%f\t", buffer[i]);
  }
  printf("\n");
}

