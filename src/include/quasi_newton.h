/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        quasi_newton.h
 * Version:     0.1.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 08-Jul-2012.
 */

#ifndef OPTIMIZATION_QUASI_NEWTON_H
#define OPTIMIZATION_QUASI_NEWTON_H

#include "non_linear_component.h"
#include "line_search_component.h"

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

typedef struct _QuasiNewtonParameter {
    char formula;
    double tolerance;
    int upper_iter;
} QuasiNewtonParameter;

int
quasi_newton(
    double *x,
    double **b,
    int n,
    FunctionObject *function_object,
    line_search_t line_search,
    LineSearchParameter *line_search_parameter,
    QuasiNewtonParameter *quasi_newton_parameter
);

#endif // OPTIMIZATION_QUASI_NEWTON_H

