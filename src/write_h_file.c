#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]) {
  int i;
  int n;
  double h;
  if (argc != 3) {
    printf("usage: %s <h-constant> <n>\n", argv[0]);
    exit(-1);
  }

  h = atof(argv[1]);
  n = atoi(argv[2]);

  FILE* fp = fopen("H.in", "w");
  if (!fp) {
    printf("ERROR opening %s for writing\n", "H.in");
    exit(-1);
  }

  for (i = 0; i < n; i++) {
    fwrite(&h, 1, sizeof(double), fp);
  }

  fclose(fp);
}
