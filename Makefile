# vim:set ts=8 sts=4 sw=4 tw=0:
#
# File:        Makefile
# Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
# Last Change: 18-Jun-2012.
#
# Makefile for drivers

CC = gcc
CFLAGS = -Wall -O3

all: driver1 driver2 mv

driver1: quasi_newton_bfgs.o line_search.o non_linear_component.o mymath.o myvector.o mymatrix.o driver1.o
	$(CC) $(CFLAGS) quasi_newton_bfgs.o line_search.o non_linear_component.o mymath.o myvector.o mymatrix.o driver1.o -o driver1 -lm

driver2: quasi_newton_bfgs.o line_search.o non_linear_component.o mymath.o myvector.o mymatrix.o driver2.o
	$(CC) $(CFLAGS) quasi_newton_bfgs.o line_search.o non_linear_component.o mymath.o myvector.o mymatrix.o driver2.o -o driver2 -lm

driver1.o: driver1.c quasi_newton_bfgs.o non_linear_component.o
	$(CC) -c driver1.c

driver2.o: driver2.c quasi_newton_bfgs.o non_linear_component.o
	$(CC) -c driver2.c

quasi_newton_bfgs.o: src/quasi_newton_bfgs.c line_search.o non_linear_component.o mymath.o myvector.o mymatrix.o
	$(CC) -c src/quasi_newton_bfgs.c

line_search.o: src/line_search.c non_linear_component.o mymath.o
	$(CC) -c src/line_search.c

non_linear_component.o: src/non_linear_component.c
	$(CC) -c $^

mymath.o: src/mymath.c
	$(CC) -c $^

myvector.o: src/myvector.c
	$(CC) -c $^

mymatrix.o: src/mymatrix.c
	$(CC) -c $^

mv:
	mv *.o bin

clean:
	rm -vf bin/*.o driver1 driver2

