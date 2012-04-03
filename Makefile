EXECS = part2 driver

all: $(EXECS)

part2: part2.c
	gcc -o part2 part2.c

driver: driver.c my_dgemm.h
	gcc -lcblas -lblas -o driver driver.c

clean:
	rm $(EXECS)