/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        non_linear_component.h
 * Version:     0.1.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 22-Jun-2012.
 */

#ifndef OPTIMIZATION_NON_LINEAR_COMPONENT_H
#define OPTIMIZATION_NON_LINEAR_COMPONENT_H

enum NonLinearStatus {
    NON_LINEAR_FUNCTION_NAN = -1,
    NON_LINEAR_SATISFIED = 0,
};

typedef struct _FunctionObject {
    double  (*function)(const double *, int);
    void    (*gradient)(double *, const double *, int);
} FunctionObject;

typedef struct _NonLinearComponent {
    int iteration_f;
    int iteration_g;
    double f;
    double alpha;
    FunctionObject *function_object;
} NonLinearComponent;

void
initialize_non_linear_component(
    FunctionObject *function_object,
    NonLinearComponent *component
);

int evaluate_function(
    const double *x,
    int n,
    NonLinearComponent *component
);

int evaluate_gradient(
    double *g,
    const double *x,
    int n,
    NonLinearComponent *component
);

int evaluate_function_gradient(
    double *g,
    const double *x,
    int n,
    NonLinearComponent *component
);

#endif // OPTIMIZATION_NON_LINEAR_COMPONENT_H

