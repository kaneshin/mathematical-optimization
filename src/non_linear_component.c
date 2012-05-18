// vim:set ts=8 sts=4 sw=4 tw=0:
// vim:set foldmethod=marker foldmarker={{{,}}}:
/* ===========================================================================
 *  File: non_linear_component.c
 *  Version: 0.9.0
 *  Last Change: 18-May-2012.
 *  Maintainer: Shintaro Kaneko <kaneshin0120@gmail.com>
 *  Description:
=========================================================================== */

#include "include/non_linear_component.h"

int
evaluate_function
(
    const double *x,
    int n,
    NonLinearComponent *component
)
{
    component->f = component->function_object->function(x, n);
    component->iteration_f++;
    if (component->f != component->f) {
        return NON_LINEAR_FUNCTION_NAN;
    }
    return NON_LINEAR_SATISFIED;
}

int
evaluate_gradient
(
    double *g,
    const double *x,
    int n,
    NonLinearComponent *component
)
{
    int i;
    component->function_object->gradient(g, x, n);
    component->iteration_g++;
    for (i = 0; i < n; i++) {
        if ( g[i] != g[i] ) {
            return NON_LINEAR_FUNCTION_NAN;
        }
    }
    return NON_LINEAR_SATISFIED;
}

int
evaluate_function_gradient
(
    double *g,
    const double *x,
    int n,
    NonLinearComponent *component
)
{
    int i;
    component->f = component->function_object->function(x, n);
    component->iteration_f++;
    component->function_object->gradient(g, x, n);
    component->iteration_g++;
    if (component->f != component->f) {
        return NON_LINEAR_FUNCTION_NAN;
    }
    for (i = 0; i < n; i++) {
        if ( g[i] != g[i] ) {
            return NON_LINEAR_FUNCTION_NAN;
        }
    }
    return NON_LINEAR_SATISFIED;
}

