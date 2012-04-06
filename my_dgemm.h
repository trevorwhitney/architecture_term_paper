#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <cblas.h>

#define MIN(a,b) ( (a) < (b) ? (a) : (b))

double my_dgemm(int n, double *a, double *b, double *c);
double my_cblas_dgemm(int n, double *a, double *b, double *c);
double step01(int n, double *a, double *b, double *c, int t);

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

double my_cblas_dgemm(int n, double *a, double *b, double *c) {
  clock_t start, end;
  double elapsed_time;
  
  start = clock();
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
    n, n, n, 1, a, n, b, n, 0, c, n);
  end = clock();

  elapsed_time = ((double) (end-start))/CLOCKS_PER_SEC;
  return elapsed_time;
}

double step01(int n, double *a, double *b, double *c, int t) {
  int i, j ,k, ii, kk;
  clock_t start, end;
  double elapsed_time;
  register double r;

  start = clock();
  for (kk = 1; kk < n; kk += t) {
    for (ii = 1; ii < n; ii += t) {
      for (j = 0; j < n; j++) {
        for (k = kk; k < MIN((kk+t-1), n); k++) {
          r = b[k+n*j];
          for (i = ii; i < MIN((ii+t-1), n); i++) {
            c[i+n*j] = c[i+n*j] + a[i+n*k] * b[k+n*j];
          }
        }
      }
    }
  }
  end = clock();
  elapsed_time = ((double) (end-start))/CLOCKS_PER_SEC;
  return elapsed_time;
}