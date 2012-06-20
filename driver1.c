/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        driver1.c
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 20-Jun-2012.
 *
 * Problem:
 * 	minimize f(x) = (x1 - x2^2)^2 / 2 + (x2 - 2)^2 / 2
 *
 * 	gradient f(x) = [
 * 		gf(x)1 = x1 - x2^2
 * 		gf(x)2 = -2 * x2(x1 - x2^2) + x2 - 2
 * 	]
 *
 * 	x* = [ 4, 2 ]^T
 *
 */

#include <stdlib.h>

#include "src/include/quasi_newton_bfgs.h"
#include "src/include/non_linear_component.h"

static double
function(const double *x, int n);

static void
gradient(double *g, const double *x, int n);

int
main(int argc, char* argv[]) {
    int i, n;
    double *x, **b;
    FunctionObject Function;

    n = 2;
    x = (double *)malloc(sizeof(double) * n);
    b = (double **)malloc(sizeof(double *) * n);
    *b = (double *)malloc(sizeof(double) * n * n);
    for (i = 1; i < n; ++i) b[i] = b[i - 1] + n;

    for (i = 0; i < n; ++i) x[i] = 0.;
    b[0][0] = 1.;   b[0][1] = -2.;
    b[1][0] = -2.;  b[1][1] = 6.;

    Function.function = function;
    Function.gradient = gradient;

    quasi_newton_bfgs(x, b, n, &Function, NULL);

    if (NULL != x) {
        free(x);
        x = NULL;
    }
    if (NULL != b) {
        if (NULL != *b) {
            free(*b);
        }
        free(b);
        b = NULL;
    }
    return 0;
}

static double
function(const double *x, int n) {
    return (x[0] - x[1] * x[1]) * (x[0] - x[1] * x[1]) + (x[1] - 2.) * (x[1] - 2.) / 2.;
}

static void
gradient(double *g, const double *x, int n) {
    g[0] = x[0] - x[1] * x[1];
    g[1] = -2. * x[1] * (x[0] - x[1] * x[1]) + x[1] - 2.;
}

