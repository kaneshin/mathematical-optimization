/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        print_message.c
 * Version:     0.1.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 09-Jul-2012.
 * TODO:
 */

#include "include/print_message.h"

#include <stdio.h>

void
print_iteration_info(
    int iteration,
    double g_norm,
    NonLinearComponent *component
) {
    printf("\n == iteration: %7d ================================\n", iteration);
    printf("step width parameter:  \t%13.6e\n", component->alpha);
    printf("-------------------------------------------------------\n");
    printf("function value:        \t%13.6e\n", component->f);
    printf("||gf(x_k+1)||_infinity = %e\n", g_norm);
    printf("-------------------------------------------------------\n");
}

void
print_result_info(
    int status,
    int iteration,
    NonLinearComponent *component
) {
    printf("\n\n\nCompute status: %3d\n", status);
    switch (status) {
        case NON_LINEAR_SATISFIED:
            printf("Satisfied: %s Method is finished\n", component->method_name);
            break;
        case NON_LINEAR_FUNCTION_NAN:
            printf("Failed: Function value is Not a Number\n");
            break;
        case NON_LINEAR_OUT_OF_MEMORY:
            printf("Failed: Out of Memory\n");
            break;
        case NON_LINEAR_NO_FUNCTION:
            printf("Failed: Function Object is not defined\n");
            break;
        case NON_LINEAR_NO_PARAMETER:
            printf("Failed: Parameter of line search is not defined\n");
            break;
        case NON_LINEAR_FAILED:
            printf("Failed: FAILED\n");
            break;
        case NON_LINEAR_NO_CONVERGENCE:
            printf("Failed: No convergence\n");
            break;
        case NON_LINEAR_LINE_SEARCH_FAILED:
            printf("Failed: Line Search is failed\n");
            break;
        default:
            break;
    }
    if (status >= NON_LINEAR_SATISFIED) {
        printf("-------------------------------------------------------\n");
        printf("iterations:          %12d\n", iteration);
        printf("function evaluations:%12d\n", component->iteration_f);
        printf("gradient evaluations:%12d\n", component->iteration_g);
        printf("=======================================================\n");
        printf("function value:      \t%13.6e\n", component->f);
    }
}

