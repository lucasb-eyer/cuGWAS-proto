#include "io.h"

int write_general(void* buffer, FILE* fp, int size, int start) {

  if(!fp) {
    printf("error: fp null\n");
  }
  fseek(fp, start*sizeof(double), SEEK_SET);
  //just write bytes
  size_t out = fwrite(buffer, sizeof(double), size, fp);
  printf("writing:\n");
  print_buffer(buffer, size);
  if(out != size) {
    printf("error: actual write size != proposed write size");
  }
  return 0;
}
int write_float(float *buffer, FILE* file, int length, int num_indices, int start_index) {
  return write_general((void*)buffer, file, num_indices*length*sizeof(float), start_index*length*sizeof(float));
}
int write_double(double *buffer, FILE* file, int length, int num_indices, int start_index) {
  return write_general((void*)buffer, file, num_indices*length, start_index*length);
}

int read_general(void* buffer, FILE* fp, int size, int start) {

  fseek(fp, start*sizeof(double), SEEK_SET);
  //just write bytes
  size_t out = fread(buffer, sizeof(double), size, fp);
  if(out != size) {
    printf("error: actual read size != proposed read size");
  }
  return 0;
}
int read_float(float *buffer, FILE* file, int length, int num_indices, int start_index) {
  return read_general((void*)buffer, file, num_indices*length*sizeof(float), start_index*length*sizeof(float));
}
int read_double(double *buffer, FILE* file, int length, int num_indices, int start_index) {
  return read_general((void*)buffer, file, num_indices*length, start_index*length);
}

void print_buffer(double *buffer, int items) {
  int i;
  for(i = 0; i < items; i++) {
    printf("%f\t", buffer[i]);
  }
  printf("\n");
}

