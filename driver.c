#include "my_dgemm.h"

int main (int argc, int *argv[]) {
  //matrices
  double* a;
  double* b;
  double* c;

  int i, j, k, l, size, t;
  double elapsed_time;
  //FILE *my_csv, *atlas_csv;
  FILE *results;

  //my_csv = fopen("my_dgemm.csv", "a+");
  //atlas_csv = fopen("atlas.csv", "a+");
  //fprintf(my_csv, "Iteration #, Size of Array, Time (Seconds)\n");
  //fprintf(atlas_csv, "Iteration #, Size of Array, Time (Seconds)\n");

  results = fopen("tiles.csv", "a+");
  fprintf(results, "Iteration #, Size of Array, Size of Tiles, Time (Seconds)\n");


  for (i = 1; i <= 10; i++) {
    printf("Starting iteration %d\n", i);
    size = i * 500;
    a = (double*) malloc(size*size*sizeof(double));
    b = (double*) malloc(size*size*sizeof(double));
    c = (double*) malloc(size*size*sizeof(double));

    for (j = 0; j < size * size; j++) {
      a[j] = (j % 64);
      b[j] = (j % 64);
    }

    for(t = 2; t < 250; t++) {
      elapsed_time = step01(size, a, b, c, t);
      fprintf(results, "%d, %d, %d, %.2f\n", i, size, t, elapsed_time);
    }

    //elapsed_time = my_dgemm(size, a, b, c);
    //fprintf(my_csv, "%d, %d, %.2f\n", i, size, elapsed_time);

    //elapsed_time = my_cblas_dgemm(size, a, b, c);
    //fprintf(atlas_csv, "%d, %d, %.2f\n", i, size, elapsed_time);

  }

  //fclose(my_csv);
  //fclose(atlas_csv);
  fclose(results);

  return 0;
}