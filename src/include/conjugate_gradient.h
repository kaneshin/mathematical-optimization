/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        conjugate_gradient.h
 * Version:     0.1.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 08-Jul-2012.
 */

#ifndef OPTIMIZATION_CONJUGATE_GRADIENT_H
#define OPTIMIZATION_CONJUGATE_GRADIENT_H

#include "non_linear_component.h"
#include "line_search_component.h"

enum ConjugateGradientStatus {
    CONJUGATE_GRADIENT_FUNCTION_NAN = -8,
    CONJUGATE_GRADIENT_OUT_OF_MEMORY,
    CONJUGATE_GRADIENT_NO_FUNCTION,
    CONJUGATE_GRADIENT_NO_PARAMETER,
    CONJUGATE_GRADIENT_SATISFIED = 0,
    CONJUGATE_GRADIENT_FAILED,
    CONJUGATE_GRADIENT_LINE_SEARCH_FAILED,
    CONJUGATE_GRADIENT_NOT_UPDATE,
};

typedef int (*line_search_t)(
    double *,
    const double *,
    const double *,
    const double *,
    int,
    EvaluateObject *,
    LineSearchParameter *,
    NonLinearComponent *
);

typedef struct _ConjugateGradientParameter {
    double tolerance;
    int upper_iter;
} ConjugateGradientParameter;

int
conjugate_gradient(
    double *x,
    int n,
    FunctionObject *function_object,
    line_search_t line_search,
    LineSearchParameter *line_search_parameter,
    ConjugateGradientParameter *conjugate_gradient_parameter
);

#endif //  OPTIMIZATION_CONJUGATE_GRADIENT_H

