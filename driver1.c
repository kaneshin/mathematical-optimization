/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        driver1.c
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 05-Jul-2012.
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
 */

#include <stdlib.h>

#include "src/include/quasi_newton_bfgs.h"
#include "src/include/non_linear_component.h"
#include "src/include/line_search_component.h"

#define __LINE_SEARCH_METHOD 1
#if __LINE_SEARCH_METHOD == 1
    #include "src/include/armijo.h"
#elif __LINE_SEARCH_METHOD == 2
    #include "src/include/wolfe.h"
#elif __LINE_SEARCH_METHOD == 3
    #include "src/include/strong_wolfe.h"
#elif __LINE_SEARCH_METHOD == 4
    #include "src/include/backtracking_wolfe.h"
#elif __LINE_SEARCH_METHOD == 5
    #include "src/include/backtracking_strong_wolfe.h"
#endif

static double
function(const double *x, int n);

static void
gradient(double *g, const double *x, int n);

int
main(int argc, char* argv[]) {
    int i, n;
    double *x, **b;
    FunctionObject Function;
    LineSearchParameter line_search_parameter;

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
#ifdef OPTIMIZATION_LINE_SEARCH_ARMIJO_H
    default_armijo_parameter(&line_search_parameter);
#endif
#ifdef OPTIMIZATION_LINE_SEARCH_WOLFE_H
    default_wolfe_parameter(&line_search_parameter);
#endif
#ifdef OPTIMIZATION_LINE_SEARCH_STRONG_WOLFE_H
    default_strong_wolfe_parameter(&line_search_parameter);
#endif
#ifdef OPTIMIZATION_LINE_SEARCH_BACKTRACKING_WOLFE_H
    default_backtracking_wolfe_parameter(&line_search_parameter);
#endif
#ifdef OPTIMIZATION_LINE_SEARCH_BACKTRACKING_STRONG_WOLFE_H
    default_backtracking_strong_wolfe_parameter(&line_search_parameter);
#endif

    /* int
     * quasi_newton_bfgs(
     *     double *x,
     *     double **b,
     *     int n,
     *     FunctionObject *function_object,
     *     line_search_t line_search,
     *     LineSearchParameter *line_search_parameter,
     *     QuasiNewtonBFGSParameter *quasi_newton_bfgs_parameter
     * )
     */
    quasi_newton_bfgs(
            x,
            b,
            n,
            &Function,
#ifdef OPTIMIZATION_LINE_SEARCH_ARMIJO_H
            armijo
#endif
#ifdef OPTIMIZATION_LINE_SEARCH_WOLFE_H
            wolfe
#endif
#ifdef OPTIMIZATION_LINE_SEARCH_STRONG_WOLFE_H
            strong_wolfe
#endif
#ifdef OPTIMIZATION_LINE_SEARCH_BACKTRACKING_WOLFE_H
            backtracking_wolfe
#endif
#ifdef OPTIMIZATION_LINE_SEARCH_BACKTRACKING_STRONG_WOLFE_H
            backtracking_strong_wolfe
#endif
            ,
            &line_search_parameter,
            NULL
    );

    if (NULL != x) {
        free(x);
        x = NULL;
    }
    if (NULL != b) {
        if (NULL != *b) {
            free(*b);
            *b = NULL;
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

