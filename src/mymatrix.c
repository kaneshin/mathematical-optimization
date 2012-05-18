// vim:set ts=8 sts=4 sw=4 tw=0:
// vim:set foldmethod=marker foldmarker={{{,}}}:
/* ===========================================================================
 *  File: mymatrix.c
 *  Version: 0.9.0
 *  Last Change: 18-May-2012.
 *  Maintainer: Shintaro Kaneko <kaneshin0120@gmail.com>
 *  Description:
=========================================================================== */

#include "include/mymatrix.h"

#include <stdlib.h>

void
matrix_delete(double **matrix)
{
    if (NULL != matrix) {
        if (NULL != *matrix) {
            free(*matrix);
        }
        free(matrix);
        matrix = NULL;
    }
}

void
identity_matrix(double **a, int n)
{
    int i, j;
    for (i = 0; i < n; i++) {
        a[i][i] = 1.;
        for (j = 0; j < i; j++) {
            a[i][j] = 0.;
        }
        for (j = i + 1; j < n; j++) {
            a[i][j] = 0.;
        }
    }
}

