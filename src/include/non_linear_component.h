/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        non_linear_component.h
 * Version:     0.1.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 30-Jun-2012.
 */

#ifndef OPTIMIZATION_NON_LINEAR_COMPONENT_H
#define OPTIMIZATION_NON_LINEAR_COMPONENT_H

enum NonLinearStatus {
    NON_LINEAR_FUNCTION_NAN = -1,
    NON_LINEAR_SATISFIED = 0,
};

typedef struct _FunctionObject {
    double  (*function)(const double *, unsigned int);
    void    (*gradient)(double *, const double *, unsigned int);
} FunctionObject;

typedef struct _NonLinearComponent {
    unsigned int iteration_f;
    unsigned int iteration_g;
    double f;
    double alpha;
    FunctionObject *function_object;
} NonLinearComponent;

typedef struct _EvaluateObject {
    int (*function)(
            const double *,
            unsigned int,
            NonLinearComponent *
        );
    int (*gradient)(
            double *,
            const double *,
            unsigned int,
            NonLinearComponent *
        );
    int (*function_gradient)(
            double *,
            const double *,
            unsigned int,
            NonLinearComponent *
            );
    FunctionObject *function_object;
} EvaluateObject;

void
initialize_non_linear_component(
    FunctionObject *function_object,
    EvaluateObject *evaluate_object,
    NonLinearComponent *component
);

#endif // OPTIMIZATION_NON_LINEAR_COMPONENT_H

