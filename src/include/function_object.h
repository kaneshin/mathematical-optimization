#ifndef FUNCTION_OBJECT_H
#define FUNCTION_OBJECT_H

struct _NonLinearComponent;

typedef struct _FunctionObject
{
    double (*function)(double *, int);
    void (*gradient)(double *, double *, int);
} FunctionObject;

typedef struct _EvaluateFunctionObject
{
    double (*function)(double *, int, _NonLinearComponent);
    void (*gradient)(double *, double *, int, _NonLinearComponent);
} EvaluateFunctionObject;

#endif // FUNCTION_OBJECT_H

