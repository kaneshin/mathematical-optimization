/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        line_search_component.c
 * Version:     0.1.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 30-Jun-2012.
 */

#include "include/line_search_component.h"

void
default_line_search_parameter(
    LineSearchParameter *parameter
) {
    parameter->step_width   = 1.;
    parameter->xi           = 0.001;
    parameter->tau          = .5;
    parameter->sigma        = .2;
    parameter->decreasing   = .5;
    parameter->increasing   = 2.1;
}

