/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        line_search_component.h
 * Version:     0.1.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 22-Jun-2012.
 */

#ifndef OPTIMIZATION_LINE_SEARCH_COMPONENT_H
#define OPTIMIZATION_LINE_SEARCH_COMPONENT_H

enum LineSearchStatus {
    LINE_SEARCH_FUNCTION_NAN = -5,
    LINE_SEARCH_OUT_OF_MEMORY,
    LINE_SEARCH_SATISFIED = 0,
    LINE_SEARCH_FAILED,
    LINE_SEARCH_STEP_WIDTH_FAILED,
};

typedef struct _LineSearchParameter {
    double step_width;
    double xi;
    double tau;
    double sigma;
    double decreasing;
    double increasing;
} LineSearchParameter;

void
default_line_search_parameter(
    LineSearchParameter *parameter
);

#endif // OPTIMIZATION_LINE_SEARCH_COMPONENT_H

