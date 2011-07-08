#include "io.h"

int write_general(void* buffer, FILE* fp, int size, int start) {

  if(!fp) {
    printf("error: fp null\n");
  }
  fseek(fp, start*sizeof(double), SEEK_SET);
  //just write bytes
  size_t out = fwrite(buffer, sizeof(double), size, fp);
  if(out != size) {
    printf("error: actual write size != proposed write size\n");
  }
  return 0;
}

int write(double *buffer, FILE* file, int length, int num_indices, int start_index) {
  return write_general((void*)buffer, file, num_indices*length, start_index*length);
}

int read_general(void* buffer, FILE* fp, int size, int start) {

}

int read(double *buffer, FILE* fp, int length, int num_indices, int start_index) {

  fseek(fp, start_index*length*sizeof(double), SEEK_SET);

  size_t out = fread(buffer, sizeof(double), num_indices*length, fp);
  if(out != num_indices*length) {
    printf("error: actual read size(%d) != proposed read size(%d)\n", (int)out, (int)num_indices*length);

  }

  return 0; 
}

void print_buffer(double *buffer, int items) {
  int i;
  for(i = 0; i < items; i++) {
    printf("%f\t", buffer[i]);
  }
  printf("\n");
}

