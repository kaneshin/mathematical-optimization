#ifndef OPTIMIZATION_MYMATH_H
#define OPTIMIZATION_MYMATH_H

enum MyMathStatus {
    MY_MATH_FUNCTION_NAN = -5,
    MY_MATH_OUT_OF_MEMORY,
    MY_MATH_NON_FUNCTION,
    MY_MATH_SATISFIED = 0,
    MY_MATH_FAILED,
    MY_MATH_NOT_UPDATE,
};

int gauss_seidel
(
    double **a,
    double *x,
    const double *b,
    int n,
    double epsilon
);

int successive_over_relaxation
(
    double **a,
    double *x,
    const double *b,
    int n,
    double epsilon,
    double omega
);

#endif // OPTIMIZATION_MYMATH_H

