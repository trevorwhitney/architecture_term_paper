#include "my_dgemm.h"

void zero_array(int n, double *c);

int main (int argc, int *argv[]) {
  //matrices
  double* a;
  double* b;
  double* c;
  double* atlas_c;

  int i, j, k, l, size, t;
  double elapsed_time, error;
  FILE *results;

  results = fopen("results.csv", "a+");
  fprintf(results, "Step #, Size of Array, Size of Tiles, Time (Seconds), Error\n");


  for (i = 1; i <= 1; i++) {
    printf("Starting iteration %d\n", i);
    size = i * 5;
    a = (double*) malloc(size*size*sizeof(double));
    b = (double*) malloc(size*size*sizeof(double));
    c = (double*) malloc(size*size*sizeof(double));
    atlas_c = (double*) malloc(size*size*sizeof(double));

    for (j = 0; j < size * size; j++) {
      a[j] = (j % 64)*1.0;
      b[j] = (j % 64)*1.0;
    }

    for (j = 0; j < size * size; j++) {
      printf("a[%d] = %.2f\n", j, a[j]);
      printf("b[%d] = %.2f\n", j, b[j]);
    }

    //set tile size
    t = 9;

    //baseline
    elapsed_time = my_cblas_dgemm(size, a, b, atlas_c);
    error = calculate_error(atlas_c, atlas_c, size);
    fprintf(results, "atlas, %d, %d, %.2f, %.30f\n", size, t, elapsed_time, error);
    for (j = 0; j < size * size; j++) {
      printf("atlas_c[%d] = %.2f\n", j, atlas_c[j];
    }
    
    //my_dgemm
    elapsed_time = my_dgemm(size, a, b, c);
    error = calculate_error(c, atlas_c, size);
    fprintf(results, "my_dgemm, %d, %d, %.2f, %.30f\n", size, t, elapsed_time, error);
    for (j = 0; j < size * size; j++) {
      printf("my_dgemm[%d] = %.2f\n", j, c[j];
    }
    zero_array(size, c);

    //step01
    elapsed_time = step01(size, a, b, c, t);
    error = calculate_error(c, atlas_c, size);
    fprintf(results, "step01, %d, %d, %.2f, %.30f\n", size, t, elapsed_time, error);
    for (j = 0; j < size * size; j++) {
      printf("step01[%d] = %.2f\n", j, c[j];
    }
    zero_array(size, c);

    //step02
    elapsed_time = step02(size, a, b, c, t);
    error = calculate_error(c, atlas_c, size);
    fprintf(results, "step02, %d, %d, %.2f, %.30f\n", size, t, elapsed_time, error);
    for (j = 0; j < size * size; j++) {
      printf("step02[%d] = %.2f\n", j, c[j];
    }
    zero_array(size, c);

    //step03
    elapsed_time = step03(size, a, b, c, t);
    error = calculate_error(c, atlas_c, size);
    fprintf(results, "step03, %d, %d, %.2f, %.30f\n", size, t, elapsed_time, error);
    for (j = 0; j < size * size; j++) {
      printf("step03[%d] = %.2f\n", j, c[j];
    }
    zero_array(size, c);

    //step04
    elapsed_time = step04(size, a, b, c, t);
    error = calculate_error(c, atlas_c, size);
    fprintf(results, "step04, %d, %d, %.2f, %.30f\n", size, t, elapsed_time, error);
    for (j = 0; j < size * size; j++) {
      printf("step04[%d] = %.2f\n", j, c[j];
    }
  }

  fclose(results);

  return 0;
}

void zero_array(int n, double *c) {
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      c[i + n*j] = 0.0;
}