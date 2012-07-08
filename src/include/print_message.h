/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        print_message.h
 * Version:     0.1.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 08-Jul-2012.
 */

#ifndef OPTIMIZATION_PRINT_MESSAGE_H
#define OPTIMIZATION_PRINT_MESSAGE_H

#include "non_linear_component.h"

void
print_iteration_info(
    int iteration,
    double g_norm,
    NonLinearComponent *component
);

void
print_result_info(
    int status,
    int iteration,
    NonLinearComponent *component
);

#endif // OPTIMIZATION_PRINT_MESSAGE_H

