// vim:set ts=8 sts=4 sw=4 tw=0:
// vim:set foldmethod=marker foldmarker={{{,}}}:
/* ===========================================================================
 *  File: mymath.c
 *  Version: 0.9.0
 *  Last Change: 18-May-2012.
 *  Maintainer: Shintaro Kaneko <kaneshin0120@gmail.com>
 *  Description:
=========================================================================== */

#include "include/mymath.h"

#include <math.h>

int
gauss_seidel
(
    double **a,
    double *x,
    const double *b,
    int n,
    double epsilon
)
{
    int i, j;
    double norm, temp, x_old;
    do {
        norm = 0.;
        for (i = 0; i < n; i++) {
            x_old = x[i];
            temp = b[i];
            for (j = 0; j < i; j++) {
                temp -= a[i][j] * x[j];
            }
            for (j = i + 1; j < n; j++) {
                temp -= a[i][j] * x[j];
            }
            if (0 != a[i][i]) {
                x[i] = temp / a[i][i];
            }
            temp = fabs(x[i] - x_old);
            if (norm < temp) norm = temp;
        }
    } while(norm > epsilon);
    return MY_MATH_SATISFIED;
}

int
successive_over_relaxation
(
    double **a,
    double *x,
    const double *b,
    int n,
    double epsilon,
    double omega
)
{
    int i, j;
    double norm, temp, x_old;
    do {
        norm = 0.;
        for (i = 0; i < n; i++) {
            x_old = x[i];
            temp = b[i];
            for (j = 0; j < i; j++) {
                temp -= a[i][j] * x[j];
            }
            for (j = i + 1; j < n; j++) {
                temp -= a[i][j] * x[j];
            }
            if (0 != a[i][i]) {
                x[i] = x_old + omega * (temp / a[i][i] - x_old);
            }
            temp = fabs(x[i] - x_old);
            if(norm < temp) norm = temp;
        }
    } while(norm > epsilon);
    return MY_MATH_SATISFIED;
}

