#ifndef OPTIMIZATION_LINE_SEARCH_BACKTRACKING_WOLFE_H
#define OPTIMIZATION_LINE_SEARCH_BACKTRACKING_WOLFE_H

#include "non_linear_component.h"
#include "line_search.h"

int
backtracking_wolfe(
    double *x,
    double *g,
    double *d,
    int n,
    LineSearchParameter *line_search_parameter,
    NonLinearComponent *component
);

#endif // OPTIMIZATION_LINE_SEARCH_BACKTRACKING_WOLFE_H

