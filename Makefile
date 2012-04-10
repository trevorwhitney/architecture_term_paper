EXECS = part2 driver

all: $(EXECS)

part2: part2.c
	gcc -o part2 part2.c

driver: driver.c my_dgemm.h
	gcc -lcblas -lblas -o driver driver.c

clean:
	rm $(EXECS)

optimizations: driver.c my_dgemm.h
	gcc -O1 -lcblas -lblas -o driver_O1 driver.c
	gcc -O2 -lcblas -lblas -o driver_O2 driver.c
	gcc -O3 -lcblas -lblas -o driver_O3 driver.c
	gcc -O3 -funroll-loops -lcblas -lblas -o driver_O3_funroll driver.c

tiles: driver.c my_dgemm.h
	gcc -lcblas -lblas -o tiles driver.c

test_steps: test_steps.c my_dgemm.h
	gcc -lcblas -lblas -o test_steps test_steps.c