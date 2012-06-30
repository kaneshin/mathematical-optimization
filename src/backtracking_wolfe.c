/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        backtracking_wolfe.c
 * Version:     0.1.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 30-Jun-2012.
 */

#include "include/backtracking_wolfe.h"

#include <stdlib.h>
#include <math.h>

#include "include/myvector.h"

void
default_backtracking_wolfe_parameter(
    LineSearchParameter *parameter
) {
    parameter->step_width   = 1.;
    parameter->xi           = 0.001;
    parameter->tau          = .5;
    parameter->sigma        = .2;
    parameter->decreasing   = .5;
    parameter->increasing   = 2.1;
}

int
backtracking_wolfe(
    double *x,
    double *g,
    double *d,
    unsigned int n,
    LineSearchParameter *line_search_parameter,
    EvaluateObject *evaluate_object,
    NonLinearComponent *non_linear_component
) {
    unsigned int i;
    double width, beta, f_x, gd, *x_temp, *g_temp, *storage;

    /* allocate memory to storage for x_temp and g_temp */
    if (NULL == (storage = (double *)malloc(sizeof(double) * 2 * n))) {
        return LINE_SEARCH_OUT_OF_MEMORY;
    }
    x_temp = storage;
    g_temp = x_temp + n;

    width = line_search_parameter->tau;
    beta = line_search_parameter->step_width > 0. ? line_search_parameter->step_width : 1.;
    if (NON_LINEAR_FUNCTION_NAN == evaluate_object->function(x, n, non_linear_component)) {
        return LINE_SEARCH_FUNCTION_NAN;
    }
    f_x = non_linear_component->f;
    gd = dot_product(g, d, n);
    for (i = 0; i < 20000; ++i) {
        update_step_vector(x_temp, x, beta, d, n);
        if (NON_LINEAR_FUNCTION_NAN == evaluate_object->function(x_temp, n, non_linear_component)) {
            return LINE_SEARCH_FUNCTION_NAN;
        }
        if (non_linear_component->f <= f_x + line_search_parameter->xi * beta * gd) {
            evaluate_object->gradient(g_temp, x_temp, n, non_linear_component);
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

