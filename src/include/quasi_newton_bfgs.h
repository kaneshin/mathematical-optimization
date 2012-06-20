#ifndef OPTIMIZATION_QUASI_NEWTON_BFGS_H
#define OPTIMIZATION_QUASI_NEWTON_BFGS_H

#include "non_linear_component.h"
#include "line_search.h"

enum QuasiNewtonBFGSStatus {
    QUASI_NEWTON_BFGS_FUNCTION_NAN = -5,
    QUASI_NEWTON_BFGS_OUT_OF_MEMORY,
    QUASI_NEWTON_BFGS_NON_FUNCTION,
    QUASI_NEWTON_BFGS_SATISFIED = 0,
    QUASI_NEWTON_BFGS_FAILED,
    QUASI_NEWTON_BFGS_LINE_SEARCH_FAILED,
    QUASI_NEWTON_BFGS_NOT_UPDATE,
};

typedef struct _QuasiNewtonBFGSParameter
{
    double tolerance;
    int max_iteration;
} QuasiNewtonBFGSParameter;

int quasi_newton_bfgs
(
    double *x,
    double **b,
    int n,
    FunctionObject *function_object,
    QuasiNewtonBFGSParameter *quasi_newton_bfgs_parameter,
    int (*line_search)(double *, double *, double *, int,
        LineSearchParameter *line_search_parameter,
        NonLinearComponent *component
    ),
    LineSearchParameter *line_search_parameter
);

#endif  // OPTIMIZATION_QUASI_NEWTON_METHOD_BFGS_H

