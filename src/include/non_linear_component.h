#ifndef OPTIMIZATION_NON_LINEAR_COMPONENT_H
#define OPTIMIZATION_NON_LINEAR_COMPONENT_H

typedef unsigned int uint;

enum NonLinearStatus {
    NON_LINEAR_FUNCTION_NAN = -1,
    NON_LINEAR_SATISFIED = 0,
};

typedef struct _FunctionObject
{
    double  (*function)(const double *, int);
    void    (*gradient)(double *, const double *, int);
} FunctionObject;

typedef struct _NonLinearComponent
{
    int iteration_f;
    int iteration_g;
    double f;
    double alpha;
    struct _FunctionObject *function_object;
} NonLinearComponent;

int evaluate_function
(
    const double *x,
    int n,
    NonLinearComponent *component
);

int evaluate_gradient
(
    double *g,
    const double *x,
    int n,
    NonLinearComponent *component
);

int evaluate_function_gradient
(
    double *g,
    const double *x,
    int n,
    NonLinearComponent *component
);

#endif OPTIMIZATION_NON_LINEAR_COMPONENT_H

