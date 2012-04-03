#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <cblas.h>

double my_dgemm(int n, double *a, double *b, double *c);
double my_cblas_dgemm(int n, double *a, double *b, float *c);

double my_dgemm (int n, double *a, double *b, double *c) {
  //iterators
  int i, j ,k;
  clock_t start, end;
  double elapsed_time;

  start = clock();
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      for (k = 0; k < n; k++) {
        c[i+n*j] = c[i+n*j] + a[i+n*k] * b[k+n*j];
      }
    }
  }
  end = clock();
  elapsed_time = ((double) (end-start))/CLOCKS_PER_SEC;
  return elapsed_time;
}

double my_cblas_dgemm(int n, double *a, double *b, float *c) {
  clock_t start, end;
  double elapsed_time;
  
  start = clock();
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
    n, n, n, 1, a, n, b, n, 0, *c, n);
  end = clock();

  elapsed_time = ((double) (end-start))/CLOCKS_PER_SEC;
  return elapsed_time;
}