/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        quasi_newton_bfgs.h
 * Version:     0.1.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 30-Jun-2012.
 */

#ifndef OPTIMIZATION_QUASI_NEWTON_BFGS_H
#define OPTIMIZATION_QUASI_NEWTON_BFGS_H

#include "non_linear_component.h"
#include "line_search_component.h"

enum QuasiNewtonBFGSStatus {
    QUASI_NEWTON_BFGS_FUNCTION_NAN = -8,
    QUASI_NEWTON_BFGS_OUT_OF_MEMORY,
    QUASI_NEWTON_BFGS_NO_FUNCTION,
    QUASI_NEWTON_BFGS_NO_PARAMETER,
    QUASI_NEWTON_BFGS_SATISFIED = 0,
    QUASI_NEWTON_BFGS_FAILED,
    QUASI_NEWTON_BFGS_LINE_SEARCH_FAILED,
    QUASI_NEWTON_BFGS_NOT_UPDATE,
};

typedef int (*line_search_t)(
    double *,
    double *,
    double *,
    double *,
    unsigned int,
    LineSearchParameter *,
    EvaluateObject *,
    NonLinearComponent *
);

typedef struct _QuasiNewtonBFGSParameter {
    char formula;
    double tolerance;
    int upper_iteration;
} QuasiNewtonBFGSParameter;

int
quasi_newton_bfgs(
    double *x,
    double **b,
    unsigned int n,
    FunctionObject *function_object,
    line_search_t line_search,
    LineSearchParameter *line_search_parameter,
    QuasiNewtonBFGSParameter *quasi_newton_bfgs_parameter
);

#endif  // OPTIMIZATION_QUASI_NEWTON_METHOD_BFGS_H

