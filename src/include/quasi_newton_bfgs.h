/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        quasi_newton_bfgs.h
 * Version:     0.1.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 22-Jun-2012.
 */

#ifndef OPTIMIZATION_QUASI_NEWTON_BFGS_H
#define OPTIMIZATION_QUASI_NEWTON_BFGS_H

#include "non_linear_component.h"
#include "line_search_component.h"

enum QuasiNewtonBFGSStatus {
    QUASI_NEWTON_BFGS_FUNCTION_NAN = -8,
    QUASI_NEWTON_BFGS_OUT_OF_MEMORY,
    QUASI_NEWTON_BFGS_NO_FUNCTION,
    QUASI_NEWTON_BFGS_SATISFIED = 0,
    QUASI_NEWTON_BFGS_FAILED,
    QUASI_NEWTON_BFGS_LINE_SEARCH_FAILED,
    QUASI_NEWTON_BFGS_NOT_UPDATE,
};

typedef struct _QuasiNewtonBFGSParameter {
    double tolerance;
    int upper_iteration;
} QuasiNewtonBFGSParameter;

typedef struct _BFGSFormula {
    int (*direction_search)(double *, double **, double *, int);
    int (*update_bfgs)(double **, const double *, const double *, double *, int);
} BFGSFormula;

int quasi_newton_bfgs(
    double *x,
    double **b,
    int n,
    FunctionObject *function_object,
    QuasiNewtonBFGSParameter *quasi_newton_bfgs_parameter,
    char formula,
    LineSearchParameter *line_search_parameter,
    int (*line_search)(double *, double *, double *, int,
        LineSearchParameter *, NonLinearComponent *
    )
);

#endif  // OPTIMIZATION_QUASI_NEWTON_METHOD_BFGS_H

