/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        backtracking_wolfe.h
 * Version:     0.1.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 30-Jun-2012.
 */

#ifndef OPTIMIZATION_LINE_SEARCH_BACKTRACKING_WOLFE_H
#define OPTIMIZATION_LINE_SEARCH_BACKTRACKING_WOLFE_H

#include "non_linear_component.h"
#include "line_search_component.h"

void
default_backtracking_wolfe_parameter(
    LineSearchParameter *parameter
);

int
backtracking_wolfe(
    double *x,
    double *g,
    double *d,
    unsigned int n,
    LineSearchParameter *line_search_parameter,
    EvaluateObject *evaluate_object,
    NonLinearComponent *non_linear_component
);

#endif // OPTIMIZATION_LINE_SEARCH_BACKTRACKING_WOLFE_H

