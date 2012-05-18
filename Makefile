#===========================================================================
# File: Makefile
# Last Change: 17-May-2012.
# Maintainer: Shintaro Kaneko <kaneshin0120@gmail.com>
#===========================================================================
#
# Makefile for function1.c
#

CC = gcc
CFLAGS = -Wall -O3
OBJ = bin/driver.out

driver:
	$(CC) driver.c src/non_linear_component.c src/quasi_newton_bfgs.c src/line_search.c src/mymath.c src/myvector.c src/mymatrix.c -o $(OBJ)

clean:
	rm -f $(OBJ)

