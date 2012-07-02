# vim:set ts=8 sts=4 sw=4 tw=0:
#
# File:         Makefile
# Maintainer:   Shintaro Kaneko <kaneshin0120@gmail.com>
# Last Change:  02-Jul-2012.
#
# Makefile for drivers

CC = gcc
CFLAGS = -Wall -O3
_SRCS = quasi_newton_bfgs.c\
	armijo.c\
	wolfe.c\
	strong_wolfe.c\
	backtracking_wolfe.c\
	backtracking_strong_wolfe.c\
	line_search_component.c\
	non_linear_component.c\
	mymath.c
_OBJS = $(_SRCS:%.c=%.o)
SRCDIR = src
OBJDIR = bin
SRCS = $(patsubst %,$(SRCDIR)/%,$(_SRCS))
OBJS = $(patsubst %,$(OBJDIR)/%,$(_OBJS))
PROGS = driver1 driver2

all: $(PROGS)

#####	drivers
driver%: $(OBJDIR)/driver%.o $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ -lm

#####	objects
$(OBJDIR)/driver%.o: driver%.c
	$(CC) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) -c $< -o $@

#####	post-processor
clean:
	rm -vf $(OBJS) $(PROGS)

