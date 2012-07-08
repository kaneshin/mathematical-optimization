/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        conjugate_gradient.c
 * Version:     0.1.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 08-Jul-2012.
 * TODO:
 */

#include "include/conjugate_gradient.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "include/mymath.h"
#include "include/print_message.h"

static char *method_name = "Conjugate Gradient";

static void
default_conjugate_gradient_parameter(
    ConjugateGradientParameter *parameter
);

int
conjugate_gradient(
    double *x,
    int n,
    FunctionObject *function_object,
    line_search_t line_search,
    LineSearchParameter *line_search_parameter,
    ConjugateGradientParameter *conjugate_gradient_parameter
) {
    int status;
    int i, j, iter;
    long int memory_size;
    double *storage, *storage_x,
           *d, *g, *x_temp, *g_temp, g_norm, *work, beta, g_norm_temp;
    NonLinearComponent component;
    ConjugateGradientParameter _conjugate_gradient_parameter;
    EvaluateObject evaluate_object;

    memory_size = sizeof(double) * n;
    /* allocate memory to storage for x_temp and g_temp */
    if (NULL == x) {
        /* allocate memory to storage_x */
        if (NULL == (storage_x = (double *)malloc(memory_size))) {
            status = NON_LINEAR_OUT_OF_MEMORY;
            goto result;
        }
        for (i = 0; i < n; ++i) {
            storage_x[i] = 0.;
        }
        x = storage_x;
    } else {
        storage_x = NULL;
    }
    /* allocate memory to storage for d, g, x_temp, g_temp, s, y and work */
    if (NULL == (storage = (double *)malloc(memory_size * 6))) {
        status = NON_LINEAR_OUT_OF_MEMORY;
        goto result;
    }
    d = storage;
    g = d + n;
    x_temp = g + n;
    g_temp = x_temp + n;
    work = g_temp + n;

    /* make sure that f and gf of this problem exist */
    if (NULL == function_object->function || NULL == function_object->gradient) {
        status = NON_LINEAR_NO_FUNCTION;
        goto result;
    }

    /* set the component of Non-Linear Programming */
    initialize_non_linear_component(
            method_name, function_object, &evaluate_object, &component);

    /* TODO:
     * Set defaule parameter for line_search_parameter
     * Set defaule parameter in default_conjugate_gradient_parameter
     *
     * set the parameter of Conjugate Gradient method */
    if (NULL == conjugate_gradient_parameter) {
        conjugate_gradient_parameter = &_conjugate_gradient_parameter;
        default_conjugate_gradient_parameter(conjugate_gradient_parameter);
    }

    /* parameter of Line Search */
    if (NULL == line_search_parameter) {
        status = NON_LINEAR_NO_PARAMETER;
        goto result;
    }

    /* start to compute for solving this problem */
    if (NON_LINEAR_FUNCTION_OBJECT_NAN == evaluate_object.gradient(g, x, n, &component)) {
        status = NON_LINEAR_FUNCTION_NAN;
        goto result;
    }
    /* compute initial vector of direction */
    for (i = 0; i < n; ++i) {
        d[i] = -g[i];
    }
    for (iter = 1; iter <= conjugate_gradient_parameter->upper_iter; ++iter) {
        /* compute step width with a line search algorithm */
        switch (line_search(work, x, g, d, n,
                    &evaluate_object, line_search_parameter, &component)) {
            case LINE_SEARCH_FUNCTION_NAN:
                status = NON_LINEAR_FUNCTION_NAN;
                goto result;
            case LINE_SEARCH_FAILED:
                status = NON_LINEAR_LINE_SEARCH_FAILED;
                goto result;
            default:
                break;
        }
        /* update x_temp = x + alpha * d */
        for (i = 0; i < n; ++i) {
            x_temp[i] = x[i] + component.alpha * d[i];
        }
        /* update g_temp = gradient(x_temp) */
        if (NON_LINEAR_FUNCTION_OBJECT_NAN
                == evaluate_object.gradient(g_temp, x_temp, n, &component)) {
            status = NON_LINEAR_FUNCTION_NAN;
            goto result;
        }
        /* compute infinity-norm of gradient */
        g_norm = infinity_norm(g_temp, n);

        print_iteration_info(iter, g_norm, &component);

        if (g_norm < conjugate_gradient_parameter->tolerance) {
            status = NON_LINEAR_SATISFIED;
            goto result;
        }

        /* TODO: Optimization of computing of norm of g and g_temp
         * compute beta */
        g_norm = 0.;
        g_norm_temp = 0.;
        for (i = 0; i < n; ++i) {
            g_norm += g[i] * g[i];
            g_norm_temp += g_temp[i] * g_temp[i];
        }
        beta = g_norm_temp / g_norm;
        /* update direction of descent */
        for (i = 0; i < n; ++i) {
            d[i] = -g_temp[i] + beta * d[i];
        }

        memcpy(x, x_temp, memory_size);
        memcpy(g, g_temp, memory_size);

    }

result:
    print_result_info(status, iter, &component);

    /* release memory of storage_x, and storage */
    if (NULL != storage_x) {
        free(storage_x);
        storage_x = NULL;
    }
    if (NULL != storage) {
        free(storage);
        storage = NULL;
    }

    return status;
}

static void
default_conjugate_gradient_parameter(
    ConjugateGradientParameter *parameter
) {
    parameter->tolerance = 1.e-8;
    parameter->upper_iter = 5000;
}

