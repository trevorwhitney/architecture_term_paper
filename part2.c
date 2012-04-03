#include <stdlib.h>
#include <stdio.h>
#include <time.h>

double my_dgemm (int, float*, float*, float*);

int main (int argc, int *argv[]) {
	//matrices
	float* a;
	float* b;
	float* c;
	int i, j, k, l, size;
	double elapsed_time;
	FILE *csv;

	csv = fopen("results.csv", "a+");
	fprintf(csv, "Iteration #, Size of Array, Time (Seconds)\n");

	for (i = 1; i <= 10; i++) {
		printf("Starting iteration %d\n", i);
		size = i * 500;
		a = (float*) malloc(size*size*sizeof(float));
		b = (float*) malloc(size*size*sizeof(float));
		c = (float*) malloc(size*size*sizeof(float));

		for (j = 0; j < size * size; j++) {
			a[j] = (j % 64) * 0.69;
			b[j] = (j % 64) * 0.77;
		}

		elapsed_time = my_dgemm(size, a, b, c);
		printf("Iteration %d, with an array of size %d took %.2f seconds\n", i, size, elapsed_time);
		fprintf(csv, "%d, %d, %.2f\n", i, size, elapsed_time);
	}

	fclose(csv);

	return 0;
}

double my_dgemm (int n, float *a, float *b, float *c) {
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