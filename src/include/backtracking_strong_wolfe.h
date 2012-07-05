/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        backtracking_strong_wolfe.h
 * Version:     0.2.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 02-Jul-2012.
 */

#ifndef OPTIMIZATION_LINE_SEARCH_BACKTRACKING_STRONG_WOLFE_H
#define OPTIMIZATION_LINE_SEARCH_BACKTRACKING_STRONG_WOLFE_H

#include "non_linear_component.h"
#include "line_search_component.h"

void
default_backtracking_strong_wolfe_parameter(
    LineSearchParameter *parameter
);

int
backtracking_strong_wolfe(
    double *storage,
    const double *x,
    const double *g,
    const double *d,
    int n,
    EvaluateObject *evaluate_object,
    LineSearchParameter *line_search_parameter,
    NonLinearComponent *component
);

#endif // OPTIMIZATION_LINE_SEARCH_BACKTRACKING_STRONG_WOLFE_H

