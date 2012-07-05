/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        driver6.c
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 05-Jul-2012.
 *
 * Problem:     Rastrigin function
 */

#define PI 3.14159265358979323846264338327

#include <stdlib.h>
#include <math.h>

#include "src/include/quasi_newton_bfgs.h"
#include "src/include/line_search_component.h"
#include "src/include/non_linear_component.h"

#define __LINE_SEARCH_METHOD 4
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
function(const double *x, unsigned int n);

static void
gradient(double *g, const double *x, unsigned int n);

int
main(int argc, char* argv[]) {
    unsigned int i, n;
    double *x;
    FunctionObject Function;
    LineSearchParameter line_search_parameter;

    n = 10;
    x = (double *)malloc(sizeof(double) * n);

    for (i = 0; i < n; ++i) x[i] = 1.;

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
            NULL,
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

    return 0;
}

static double
function(const double *x, unsigned int n) {
    unsigned int i;
    double f = 0., dot;
    for (i = 1, dot = x[0] * x[0]; i < n; ++i) {
        dot += x[i] * x[i];
    }
    for (i = 0; i < n; ++i) {
        f += 1 - cos(2 * PI * x[i]);
    }
    f = dot + 10.0 * f;
    return f;
}

static void
gradient(double *g, const double *x, unsigned int n) {
    unsigned int i;
    double temp;
    for (i = 0; i < n; ++i) {
        g[i] = 2 * x[i] + 20 * PI * sin(2 * PI * x[i]);
    }
}

