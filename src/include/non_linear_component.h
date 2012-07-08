/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        non_linear_component.h
 * Version:     0.1.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 08-Jul-2012.
 */

#ifndef OPTIMIZATION_NON_LINEAR_COMPONENT_H
#define OPTIMIZATION_NON_LINEAR_COMPONENT_H

enum NonLinearFunctionStatus {
    NON_LINEAR_FUNCTION_OBJECT_NAN = -1,
    NON_LINEAR_FUNCTION_OBJECT_SATISFIED = 0,
};

enum NonLinearStatus {
    NON_LINEAR_FUNCTION_NAN = -8,
    NON_LINEAR_OUT_OF_MEMORY,
    NON_LINEAR_NO_FUNCTION,
    NON_LINEAR_NO_PARAMETER,
    NON_LINEAR_SATISFIED = 0,
    NON_LINEAR_FAILED,
    NON_LINEAR_LINE_SEARCH_FAILED,
    NON_LINEAR_NOT_UPDATE,
};

typedef struct _FunctionObject {
    double  (*function)(const double *, int);
    void    (*gradient)(double *, const double *, int);
} FunctionObject;

typedef struct _NonLinearComponent {
    char *method_name;
    int iteration_f;
    int iteration_g;
    double f;
    double alpha;
    FunctionObject *function_object;
} NonLinearComponent;

typedef struct _EvaluateObject {
    int (*function)(
            const double *,
            int,
            NonLinearComponent *
        );
    int (*gradient)(
            double *,
            const double *,
            int,
            NonLinearComponent *
        );
    int (*function_gradient)(
            double *,
            const double *,
            int,
            NonLinearComponent *
            );
    FunctionObject *function_object;
} EvaluateObject;

void
initialize_non_linear_component(
    char *method_name,
    FunctionObject *function_object,
    EvaluateObject *evaluate_object,
    NonLinearComponent *component
);

#endif // OPTIMIZATION_NON_LINEAR_COMPONENT_H

