#===========================================================================
# File: Makefile
# Last Change: 05-Jun-2012.
# Maintainer: Shintaro Kaneko <kaneshin0120@gmail.com>
#===========================================================================
#
# Makefile for function1.c
#

CFLAGS = -Wall -O3
LFLAGS = -lm
SRCDIR = src
BINDIR = bin
OBJS = \
	$(BINDIR)/quasi_newton_bfgs.o \
	$(BINDIR)/line_search.o \
	$(BINDIR)/non_linear_component.o \
	$(BINDIR)/mymath.o \
	$(BINDIR)/mymatrix.o \
	$(BINDIR)/myvector.o
CC = gcc $(CFLAGS)

all: driver driver1

driver: driver.c $(OBJS)
	$(CC) -o $@.o $^ $(LFLAGS)

driver1: driver1.c $(OBJS)
	$(CC) -o $@.o $^ $(LFLAGS)

$(OBJS): mymatrix myvector mymath line_search non_linear_component quasi_newton_bfgs

quasi_newton_bfgs: non_linear_component line_search mymath myvector mymatrix
	$(CC) -c $(SRCDIR)/$@.c -o $(BINDIR)/$@.o

line_search: non_linear_component mymath
	$(CC) -c $(SRCDIR)/$@.c -o $(BINDIR)/$@.o

non_linear_component:
	$(CC) -c $(SRCDIR)/$@.c -o $(BINDIR)/$@.o

mymath:
	$(CC) -c $(SRCDIR)/$@.c -o $(BINDIR)/$@.o

mymatrix:
	$(CC) -c $(SRCDIR)/$@.c -o $(BINDIR)/$@.o

myvector:
	$(CC) -c $(SRCDIR)/$@.c -o $(BINDIR)/$@.o

clean:
	rm -f $(BINDIR)/*.o *.o

