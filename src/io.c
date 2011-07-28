#include "io.h"
#include <stdio.h>

void write(double* buffer, FILE* fp, size_t size, size_t start) {
  //  printf("Writing %d bytes at %d into %x from %x\n", size*sizeof(double), start*sizeof(double), fp, buffer);  
  fseek(fp, start*sizeof(double), SEEK_SET);
  size_t out = fwrite(buffer, sizeof(double), size, fp);
  if(out != size) {
    printf("error: actual write size != proposed write size\n");
  }
  return; 
}


void read(double *buffer, FILE* fp, size_t size, size_t start) {
  printf("Reading %lld bytes at %lld from %x into %x\n", size*sizeof(double), start*sizeof(double), fp, buffer);  
  fseek(fp, start*sizeof(double), SEEK_SET);
  size_t out = fread(buffer, sizeof(double), size, fp);
  if(out != size) {
    printf("error: actual read size(%d) != proposed read size(%d)\n", (int)out, (int)size);
  }
  return; 
}

void print_buffer(double *buffer, int items) {
  int i;
  for(i = 0; i < items; i++) {
    printf("%f\t", buffer[i]);
  }
  printf("\n");
}

