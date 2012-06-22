/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        backtracking_wolfe.c
 * Version:     0.1.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 22-Jun-2012.
 */

#include "include/backtracking_wolfe.h"

#include <stdlib.h>
#include <math.h>

#include "include/myvector.h"

int
backtracking_wolfe(
    double *x,
    double *g,
    double *d,
    int n,
    LineSearchParameter *line_search_parameter,
    NonLinearComponent *non_linear_component
) {
    int i;
    double width, beta, f_x, gd, *x_temp, *g_temp, *storage;

    storage = (double *)malloc(sizeof(double) * 2 * n);
    x_temp = storage;
    g_temp = x_temp + n;

    width = line_search_parameter->tau;
    beta = line_search_parameter->step_width > 0. ? line_search_parameter->step_width : 1.;
    if (NON_LINEAR_FUNCTION_NAN == evaluate_function(x, n, non_linear_component)) {
        return LINE_SEARCH_FUNCTION_NAN;
    }
    f_x = non_linear_component->f;
    gd = dot_product(g, d, n);
    for (i = 0; i < 20000; ++i) {
        update_step_vector(x_temp, x, beta, d, n);
        if (NON_LINEAR_FUNCTION_NAN == evaluate_function(x_temp, n, non_linear_component)) {
            return LINE_SEARCH_FUNCTION_NAN;
        }
        if (non_linear_component->f <= f_x + line_search_parameter->xi * beta * gd) {
            evaluate_gradient(g_temp, x_temp, n, non_linear_component);
            if (line_search_parameter->sigma * gd <= dot_product(g_temp, d, n)) {
                non_linear_component->alpha = beta;
                if (NULL != storage) {
                    free(storage);
                    storage = NULL;
                }
                return LINE_SEARCH_SATISFIED;
            } else {
                width = line_search_parameter->increasing;
            }
        } else {
            width = line_search_parameter->decreasing;
        }
        beta *= width;
    }
    if (NULL != storage) {
        free(storage);
        storage = NULL;
    }
    return LINE_SEARCH_FAILED;
}

