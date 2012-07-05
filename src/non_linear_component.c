/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        non_linear_component.c
 * Version:     0.1.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 02-Jul-2012.
 */

#include "include/non_linear_component.h"

static int
function(
    const double *x,
    int n,
    NonLinearComponent *component
);

static int
gradient(
    double *g,
    const double *x,
    int n,
    NonLinearComponent *component
);

static int
function_gradient(
    double *g,
    const double *x,
    int n,
    NonLinearComponent *component
);

void
initialize_non_linear_component(
    FunctionObject *function_object,
    EvaluateObject *evaluate_object,
    NonLinearComponent *component
) {
    component->iteration_f      = 0;
    component->iteration_g      = 0;
    component->f                = 0.;
    component->alpha            = 0.;
    component->function_object  = function_object;
    evaluate_object->function = function;
    evaluate_object->gradient = gradient;
    evaluate_object->function_gradient = function_gradient;
}

static int
function(
    const double *x,
    int n,
    NonLinearComponent *component
) {
    component->f = component->function_object->function(x, n);
    component->iteration_f++;
    if (component->f != component->f) {
        return NON_LINEAR_FUNCTION_NAN;
    }
    return NON_LINEAR_SATISFIED;
}

static int
gradient(
    double *g,
    const double *x,
    int n,
    NonLinearComponent *component
) {
    int i;
    component->function_object->gradient(g, x, n);
    component->iteration_g++;
    for (i = 0; i < n; ++i) {
        if (g[i] != g[i]) {
            return NON_LINEAR_FUNCTION_NAN;
        }
    }
    return NON_LINEAR_SATISFIED;
}

static int
function_gradient(
    double *g,
    const double *x,
    int n,
    NonLinearComponent *component
) {
    int i;
    component->f = component->function_object->function(x, n);
    component->iteration_f++;
    component->function_object->gradient(g, x, n);
    component->iteration_g++;
    if (component->f != component->f) {
        return NON_LINEAR_FUNCTION_NAN;
    }
    for (i = 0; i < n; ++i) {
        if (g[i] != g[i]) {
            return NON_LINEAR_FUNCTION_NAN;
        }
    }
    return NON_LINEAR_SATISFIED;
}

