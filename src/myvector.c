// vim:set ts=8 sts=4 sw=4 tw=0:
// vim:set foldmethod=marker foldmarker={{{,}}}:
/* ===========================================================================
 *  File: myvector.c
 *  Version: 0.9.0
 *  Last Change: 18-May-2012.
 *  Maintainer: Shintaro Kaneko <kaneshin0120@gmail.com>
 *  Description:
=========================================================================== */

#include "include/myvector.h"

#include <stdlib.h>
#include <math.h>

void
vector_delete(double *vector)
{
    if (NULL != vector) {
        free(vector);
        vector = NULL;
    }
}

void
copy_vector(double *x, const double *y, int n)
{
    int i;
    for (i = 0; i < n; i++) x[i] = y[i];
}

double
dot_product(const double *x, const double *y, int n)
{
    int i;
    double dot;
    dot = x[0] * y[0];
    for (i = 1; i < n; i++) dot += x[i] * y[i];
    return dot;
}

void
zero_vector(double *x, int n)
{
    int i;
    for (i = 0; i < n; i++) x[i] = 0.;
}

void
update_step_vector
(
    double *x_temp,
    const double *x, double alpha, const double *y, int n
)
{
    int i;
    for (i = 0; i < n; i++) x_temp[i] = x[i] + alpha * y[i];
}

double
norm_infty(const double *x, int n)
{
    int i;
    double norm, fabs_x;
    norm = fabs(x[0]);
    for (i = 1; i < n; i++) {
        fabs_x = fabs(x[i]);
        if (fabs_x > norm) norm = fabs_x;
    }
    return norm;
}

