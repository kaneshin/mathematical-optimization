/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        mymath.h
 * Version:     0.2.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 02-Jul-2012.
 */

#ifndef OPTIMIZATION_MYMATH_H
#define OPTIMIZATION_MYMATH_H

enum MyMathStatus {
    MY_MATH_FUNCTION_NAN = -5,
    MY_MATH_OUT_OF_MEMORY,
    MY_MATH_NON_FUNCTION,
    MY_MATH_SATISFIED = 0,
    MY_MATH_FAILED,
    MY_MATH_NOT_UPDATE,
};

double
dot_product(
    const double *x,
    const double *y,
    unsigned int n
);

void
update_step_vector(
    double *x_temp,
    const double *x,
    double alpha,
    const double *y,
    unsigned int n
);

double
manhattan_norm(
    const double *x,
    unsigned int n
);

double
euclidean_norm(
    const double *x,
    unsigned int n
);

double
infinity_norm(
    const double *x,
    unsigned int n
);

int gauss_seidel(
    double **a,
    double *x,
    const double *b,
    unsigned int n,
    double epsilon
);

int successive_over_relaxation(
    double **a,
    double *x,
    const double *b,
    unsigned int n,
    double epsilon,
    double omega
);

#endif // OPTIMIZATION_MYMATH_H

