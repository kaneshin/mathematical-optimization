/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        armijo.c
 * Version:     0.2.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 02-Jul-2012.
 */

#include "include/armijo.h"

#include "include/mymath.h"

void
default_armijo_parameter(
    LineSearchParameter *parameter
) {
    parameter->upper_iter = 5000;
    parameter->initial_step = .5;
    parameter->step_width = 1.;
    parameter->xi = 0.001;
    parameter->sigma = .2;
}

int
armijo(
    double *storage,
    const double *x,
    const double *g,
    const double *d,
    unsigned int n,
    EvaluateObject *evaluate_object,
    LineSearchParameter *parameter,
    NonLinearComponent *component
) {
    unsigned int iter;
    double width, beta, f_x, gd, *x_temp, *g_temp;

    x_temp = storage;
    g_temp = x_temp + n;

    width = parameter->initial_step;
    beta = parameter->step_width;
    if (NON_LINEAR_FUNCTION_NAN == evaluate_object->function(x, n, component))
        return LINE_SEARCH_FUNCTION_NAN;
    f_x = component->f;
    gd = dot_product(g, d, n);
    for (iter = 1; iter <= parameter->upper_iter; ++iter) {
        update_step_vector(x_temp, x, beta, d, n);
        if (NON_LINEAR_FUNCTION_NAN
                == evaluate_object->function(x_temp, n, component))
            return LINE_SEARCH_FUNCTION_NAN;
        if (component->f <= f_x + parameter->xi * beta * gd) {
            component->alpha = beta;
            return LINE_SEARCH_SATISFIED;
        }
        beta *= width;
    }
    return LINE_SEARCH_FAILED;
}
