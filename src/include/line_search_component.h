#ifndef OPTIMIZATION_LINE_SEARCH_H
#define OPTIMIZATION_LINE_SEARCH_H

#include "non_linear_component.h"

enum LineSearchStatus {
    LINE_SEARCH_FUNCTION_NAN = -1,
    LINE_SEARCH_SATISFIED = 0,
    LINE_SEARCH_FAILED,
    LINE_SEARCH_STEP_WIDTH_FAILED,
};

int line_search_armijo
(
    const double *x,
    const double *g,
    const double *d,
    int n,
    double step_width,
    double tau,
    double xi,
    NonLinearComponent *component
);

int
line_search_wolfe
(
    const double *x,
    const double *g,
    const double *d,
    int n,
    double step_width,
    double tau,
    double xi,
    double sigma,
    NonLinearComponent *component
);

int
line_search_strong_wolfe
(
    const double *x,
    const double *g,
    const double *d,
    int n,
    double step_width,
    double tau,
    double xi,
    double sigma,
    NonLinearComponent *component
);

int
line_search_backtracking_wolfe
(
    const double *x,
    const double *g,
    const double *d,
    int n,
    double step_width,
    double tau,
    double xi,
    double sigma,
    double decreasing,
    double increasing,
    NonLinearComponent *component
);

int
line_search_backtracking_strong_wolfe
(
    const double *x,
    const double *g,
    const double *d,
    int n,
    double step_width,
    double tau,
    double xi,
    double sigma,
    double decreasing,
    double increasing,
    NonLinearComponent *component
);

#endif // OPTIMIZATION_LINE_SEARCH_H

