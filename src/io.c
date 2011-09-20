#include "io.h"
#include <stdio.h>
#include <error.h>
#include <errno.h>

void write(double* buffer, FILE* fp, size_t size, size_t start) {
  fseek(fp, start*sizeof(double), SEEK_SET);
  printf("Writing: %f\n", buffer[0]);
  size_t out = fwrite(buffer, sizeof(double), size, fp);
  if(out != size) {
    printf("error: actual write size(%d) != proposed write size(%d)\n", (int)out, (int)size);
  }
  return; 
}


void read(double *buffer, FILE* fp, size_t size, size_t start) {
	/*printf("Reading %d bytes at %d (%p)\n", (int)size*sizeof(double), (int)start, buffer);*/
  fseek(fp, start*sizeof(double), SEEK_SET);
  size_t out = fread(buffer, sizeof(double), size, fp);
  if(out != size) {
	printf("ERROR: %s\n", strerror(errno));
    printf("start: %d\n", (int)start);
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

