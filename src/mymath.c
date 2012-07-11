/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        mymath.c
 * Version:     0.1.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 12-Jul-2012.
 */

#include "include/mymath.h"

#include <math.h>

/*
 * Libraries of Vector
 *  - dot_product
 *  - update_step_vector
 *  - manhattan_norm
 *  - euclidean_norm
 *  - infinity_norm
 */
double
dot_product(
    const double *x,
    const double *y,
    int n
) {
    int i;
    double dot;
    for (i = 1, dot = x[0] * y[0]; i < n; ++i)
        dot += x[i] * y[i];
    return dot;
}

void
update_step_vector(
    double *x_temp,
    const double *x,
    double alpha,
    const double *y,
    int n
) {
    int i;
    for (i = 0; i < n; ++i)
        x_temp[i] = x[i] + alpha * y[i];
}

double
manhattan_norm(
    const double *x,
    int n
) {
    /*
     * norm = |x_0| + |x_1| + ... + |x_n|
     */
    int i;
    double norm, fabs_x;
    for (i = 1, norm = fabs(x[0]); i < n; ++i)
        norm += fabs(x[i]);
    return norm;
}

double
euclidean_norm(
    const double *x,
    int n
) {
    /*
     * norm = sqrt(x_0 * x_0 + x_1 * x_1 + ... + x_n * x_n)
     */
    int i;
    double norm;
    for (i = 1, norm = x[0] * x[0]; i < n; ++i)
        norm += x[i] * x[i];
    return sqrt(norm);
}

double
infinity_norm(
    const double *x,
    int n
) {
    /*
     * norm = max(|x_0|, |x_1|, ..., |x_n|)
     */
    int i;
    double norm, fabs_x;
    for (i = 1, norm = fabs(x[0]); i < n; ++i) {
        fabs_x = fabs(x[i]);
        if (fabs_x > norm)
            norm = fabs_x;
    }
    return norm;
}

/*
 * Libraries of Methematical Analysis
 *  - gauss_seidel
 *  - successive_over_relaxation
 */
int
gauss_seidel(
    double **a,
    double *x,
    const double *b,
    int n,
    double epsilon
) {
    /*
     * NOTE:
     *  If you choose the value of Omega as Zero with SOR method, SOR method
     *  become Gause-Seidel method.
     */
    return successive_over_relaxation(a, x, b, n, epsilon, 1.);
}

int
successive_over_relaxation(
    double **a,
    double *x,
    const double *b,
    int n,
    double epsilon,
    double omega
) {
    /*
     * TODO:
     *  Error Handling if you get a value of singular as Zero.
     */
    int i, j;
    double norm, temp, x_old;
    do {
        for (i = 0, norm = 0.; i < n; ++i) {
            x_old = x[i];
            temp = b[i];
            for (j = 0; j < i; ++j)
                temp -= a[i][j] * x[j];
            for (j = i + 1; j < n; ++j)
                temp -= a[i][j] * x[j];
            if (0. != a[i][i])
                x[i] = x_old + omega * (temp / a[i][i] - x_old);
            else
                return MY_MATH_FAILED;
            temp = fabs(x[i] - x_old);
            if(norm < temp)
                norm = temp;
        }
    } while(norm > epsilon);
    return MY_MATH_SATISFIED;
}

