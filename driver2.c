/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        driver2.c
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 30-Jun-2012.
 *
 * Problem:
 * 	minimize f(x) =
 * 		sum[i=1 to n] {e^x_i - x_i * sqrt(i + 1)}
 *
 * 	gradient f(x) = [
 * 		gf(x)_i = e^x_i - sqrt(i + 1)
 * 	]
 */

#include <stdlib.h>
#include <math.h>

#include "src/include/quasi_newton_bfgs.h"
#include "src/include/line_search_component.h"
#include "src/include/backtracking_wolfe.h"
#include "src/include/non_linear_component.h"

static double
function(const double *x, unsigned int n);

static void
gradient(double *g, const double *x, unsigned int n);

int
main(int argc, char* argv[]) {
    unsigned int i, n;
    double *x;
    FunctionObject Function;
    LineSearchParameter line_search_parameter;

    n = 100;
    x = (double *)malloc(sizeof(double) * n);

    for (i = 0; i < n; ++i) x[i] = 1.;

    Function.function = function;
    Function.gradient = gradient;
    default_backtracking_wolfe_parameter(&line_search_parameter);

    quasi_newton_bfgs(x, NULL, n, &Function,
            backtracking_wolfe, &line_search_parameter, NULL);

    if (NULL != x) {
        free(x);
        x = NULL;
    }

    return 0;
}

static double
function(const double *x, unsigned int n) {
    unsigned int i;
    double f = 0.;
    for (i = 0; i < n; ++i) {
        f += exp(x[i]) - x[i] * sqrt(i + 1.);
    }
    return f;
}

static void
gradient(double *g, const double *x, unsigned int n) {
    unsigned int i;
    for (i = 0; i < n; ++i) {
        g[i] = exp(x[i]) - sqrt(i + 1.);
    }
}

