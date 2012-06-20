#ifndef OPTIMIZATION_LINE_SEARCH_H
#define OPTIMIZATION_LINE_SEARCH_H

enum LineSearchStatus {
    LINE_SEARCH_FUNCTION_NAN = -1,
    LINE_SEARCH_SATISFIED = 0,
    LINE_SEARCH_FAILED,
    LINE_SEARCH_STEP_WIDTH_FAILED,
};

typedef struct _LineSearchParameter {
    double tolerance;
    int max_iteration;
    double step_width;
    double xi;
    double tau;
    double sigma;
    double decreasing;
    double increasing;
} LineSearchParameter;

#endif // OPTIMIZATION_LINE_SEARCH_H

