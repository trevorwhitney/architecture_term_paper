#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <cblas.h>

#define MIN(a,b) ( (a) < (b) ? (a) : (b))

double calculate_error(double *my_c, double *atlas_c, int n);
double my_dgemm(int n, double *a, double *b, double *c);
double my_cblas_dgemm(int n, double *a, double *b, double *c);
double step01(int n, double *a, double *b, double *c, int t);
double step02(int n, double *a, double *b, double *c, int t);
double step03(int n, double *a, double *b, double *c, int t);
double step04(int n, double *a, double *b, double *c, int t);

double calculate_error(double *my_c, double *atlas_c, int n) {
  int i;
  double numerator;
  double denomenator;
  double error;
  double size;

  size = n*n;
  numerator = 0;
  denomenator = 0;
  for (i = 0; i < size; i++) {
    numerator += my_c[i] - atlas_c[i];
    denomenator += my_c[i];
  }

  error = (numerator*numerator)/(denomenator*denomenator);
  return error;
}

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
            c[i+n*j] = c[i+n*j] + r * a[i+n*k];
          }
        }
      }
    }
  }
  end = clock();
  elapsed_time = ((double) (end-start))/CLOCKS_PER_SEC;
  return elapsed_time;
}

double step02(int n, double *a, double *b, double *c, int t) {
  int i, j ,k, ii, jj, kk;
  clock_t start, end;
  double elapsed_time;
  register double r;

  start = clock();
  for (jj = 1; jj < n; jj += t) {
    for (kk = 1; kk < n; kk += t) {
      for (ii = 1; ii < n; ii += t) {
        for (j = jj; j < MIN((jj+t-1), n); j++) {
          for (k = kk; k < MIN((kk+t-1), n); k++) {
            r = b[k+n*j];
            for (i = ii; i < MIN((ii+t-1), n); i++) {
              c[i+n*j] = c[i+n*j] + r * a[i+n*k];
            }
          }
        }
      }
    }
  }
  end = clock();
  elapsed_time = ((double) (end-start))/CLOCKS_PER_SEC;
  return elapsed_time;
}

double step03(int n, double *a, double *b, double *c, int t) {
  int i, j ,k, ii, jj, kk;
  clock_t start, end;
  double elapsed_time;
  register double r;

  start = clock();
  for (jj = 1; jj < n; jj += t) {
    for (kk = 1; kk < n; kk += t) {
      for (ii = 1; ii < n; ii += t) {
        for (j = jj; j < MIN((jj+t-1), n); j++) {
          for (i = ii; i < MIN((ii+t-1), n); i++) {
            r = c[i+n*j];
            for (k = kk; k < MIN((kk+t-1), n); k++) {
              r += b[k+n*j] * a[i+n*k];
            }
            c[i+n*j] = r;
          }
        }
      }
    }
  }
  end = clock();
  elapsed_time = ((double) (end-start))/CLOCKS_PER_SEC;
  return elapsed_time;
}

double step04(int n, double *a, double *b, double *c, int t) {
  int i, j ,k, ii, jj, kk, dif;
  clock_t start, end;
  double elapsed_time;
  register double r;

  start = clock();
  for (jj = 1; jj < n; jj += t) {
    for (kk = 1; kk < n; kk += t) {
      for (ii = 1; ii < n; ii += t) {
        for (j = jj; j < MIN((jj+t-1), n); j++) {
          for (i = ii; i < MIN((ii+t-1), n); i++) {
            r = c[i+n*j];
            for (k = kk; k < MIN((kk+t-1), n); k+=8) {
              r += b[k+n*j] * a[i+n*k];
              r += b[k+n*j] * a[i+n*(k+1)];
              r += b[k+n*j] * a[i+n*(k+2)];
              r += b[k+n*j] * a[i+n*(k+3)];
              r += b[k+n*j] * a[i+n*(k+4)];
              r += b[k+n*j] * a[i+n*(k+5)];
              r += b[k+n*j] * a[i+n*(k+6)];
              r += b[k+n*j] * a[i+n*(k+7)];
            }
            if (n%8 > 0) {
              dif = n%8 - 1;
              while (dif > -1)
                r += b[k+n*j] * a[i+n*(n-dif--)];
            }
            c[i+n*j] = r;
          }
        }
      }
    }
  }
  end = clock();
  elapsed_time = ((double) (end-start))/CLOCKS_PER_SEC;
  return elapsed_time;
}