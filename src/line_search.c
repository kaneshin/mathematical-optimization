// vim:set ts=8 sts=4 sw=4 tw=0:
// vim:set foldmethod=marker foldmarker={{{,}}}:
/* ===========================================================================
 *  File: line_search.c
 *  Version: 0.9.0
 *  Last Change: 18-Jun-2012.
 *  Maintainer: Shintaro Kaneko <kaneshin0120@gmail.com>
 *  Description:
=========================================================================== */

#include "include/line_search.h"

#include <stdlib.h>
#include <math.h>

#include "include/mymath.h"

int
line_search_armijo
(
    const double *x,
    const double *g,
    const double *d,
    int n,
    double step_width,
    double tau,
    double xi,
    NonLinearComponent *component
)
{
    /*
     * Compute a line search routine with armijo condition
     * Return:
     */
    int i;
    double beta, f_x, gd, *x_temp, *g_temp;

    x_temp = (double *)malloc(sizeof(double) * n);
    g_temp = (double *)malloc(sizeof(double) * n);

    beta = step_width > 0. ? step_width : 1.;
    if (NON_LINEAR_FUNCTION_NAN == evaluate_function(x, n, component)) {
        return LINE_SEARCH_FUNCTION_NAN;
    }
    f_x = component->f;
    gd = dot_product(g, d, n);
    for (i = 0; i < 10000; i++) {
        update_step_vector(x_temp, x, beta, d, n);
        if (NON_LINEAR_FUNCTION_NAN == evaluate_function(x_temp, n, component)) {
            return LINE_SEARCH_FUNCTION_NAN;
        }
        if (component->f <= f_x + xi * beta * gd) {
            component->alpha = beta;
            free(x_temp);
            free(g_temp);
            return LINE_SEARCH_SATISFIED;
        }
        beta *= tau;
    }
    free(x_temp);
    free(g_temp);
    return LINE_SEARCH_FAILED;
}

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
)
{
    /*
     * Compute a line search routine wolfe armijo condition
     * Return:
     */
    int i;
    double beta, f_x, gd, *x_temp, *g_temp;

    x_temp = (double *)malloc(sizeof(double) * n);
    g_temp = (double *)malloc(sizeof(double) * n);

    beta = step_width > 0. ? step_width : 1.;
    if (NON_LINEAR_FUNCTION_NAN == evaluate_function(x, n, component)) {
        return LINE_SEARCH_FUNCTION_NAN;
    }
    f_x = component->f;
    gd = dot_product(g, d, n);
    for (i = 0; i < 10000; i++) {
        update_step_vector(x_temp, x, beta, d, n);
        if (NON_LINEAR_FUNCTION_NAN == evaluate_function(x_temp, n, component)) {
            return LINE_SEARCH_FUNCTION_NAN;
        }
        if (component->f <= f_x + xi * beta * gd) {
            evaluate_gradient(g_temp, x_temp, n, component);
            if ( sigma * gd <= dot_product(g_temp, d, n)) {
                component->alpha = beta;
                free(x_temp);
                free(g_temp);
                return LINE_SEARCH_SATISFIED;
            }
        }
        beta *= tau;
    }
    free(x_temp);
    free(g_temp);
    return LINE_SEARCH_FAILED;
}

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
)
{
    /*
     * Compute a line search routine wolfe armijo condition
     * Return:
     */
    int i;
    double beta, f_x, gd, *x_temp, *g_temp, gd_temp, sigma_gd;

    x_temp = (double *)malloc(sizeof(double) * n);
    g_temp = (double *)malloc(sizeof(double) * n);

    beta = step_width > 0. ? step_width : 1.;
    if (NON_LINEAR_FUNCTION_NAN == evaluate_function(x, n, component)) {
        return LINE_SEARCH_FUNCTION_NAN;
    }
    f_x = component->f;
    gd = dot_product(g, d, n);
    for (i = 0; i < 10000; i++) {
        update_step_vector(x_temp, x, beta, d, n);
        if (NON_LINEAR_FUNCTION_NAN == evaluate_function(x_temp, n, component)) {
            return LINE_SEARCH_FUNCTION_NAN;
        }
        if (component->f <= f_x + xi * beta * gd) {
            evaluate_gradient(g_temp, x_temp, n, component);
            gd_temp = dot_product(g_temp, d, n);
            sigma_gd = sigma * gd;
            if (sigma_gd <= gd_temp) {
                if (-sigma_gd <= fabs(gd_temp)) {
                    component->alpha = beta;
                    free(x_temp);
                    free(g_temp);
                    return LINE_SEARCH_SATISFIED;
                }
            }
        }
        beta *= tau;
    }
    free(x_temp);
    free(g_temp);
    return LINE_SEARCH_FAILED;
}

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
)
{
    /*
     * Compute a line search routine wolfe armijo condition
     * Return:
     */
    int i;
    double width, beta, f_x, gd, *x_temp, *g_temp;

    x_temp = (double *)malloc(sizeof(double) * n);
    g_temp = (double *)malloc(sizeof(double) * n);

    width = tau;
    beta = step_width > 0. ? step_width : 1.;
    if (NON_LINEAR_FUNCTION_NAN == evaluate_function(x, n, component)) {
        return LINE_SEARCH_FUNCTION_NAN;
    }
    f_x = component->f;
    gd = dot_product(g, d, n);
    for (i = 0; i < 10000; i++) {
        update_step_vector(x_temp, x, beta, d, n);
        if (NON_LINEAR_FUNCTION_NAN == evaluate_function(x_temp, n, component)) {
            return LINE_SEARCH_FUNCTION_NAN;
        }
        if (component->f <= f_x + xi * beta * gd) {
            evaluate_gradient(g_temp, x_temp, n, component);
            if ( sigma * gd <= dot_product(g_temp, d, n)) {
                component->alpha = beta;
                free(x_temp);
                free(g_temp);
                return LINE_SEARCH_SATISFIED;
            } else {
                width = increasing;
            }
        } else {
            width = decreasing;
        }
        beta *= width;
    }
    free(x_temp);
    free(g_temp);
    return LINE_SEARCH_FAILED;
}

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
)
{
    /*
     * Compute a line search routine wolfe armijo condition
     * Return:
     */
    int i;
    double width, beta, f_x, gd, *x_temp, *g_temp, gd_temp, sigma_gd;

    x_temp = (double *)malloc(sizeof(double) * n);
    g_temp = (double *)malloc(sizeof(double) * n);

    width = tau;
    beta = step_width > 0. ? step_width : 1.;
    if (NON_LINEAR_FUNCTION_NAN == evaluate_function(x, n, component)) {
        return LINE_SEARCH_FUNCTION_NAN;
    }
    f_x = component->f;
    gd = dot_product(g, d, n);
    for (i = 0; i < 10000; i++) {
        update_step_vector(x_temp, x, beta, d, n);
        if (NON_LINEAR_FUNCTION_NAN == evaluate_function(x_temp, n, component)) {
            return LINE_SEARCH_FUNCTION_NAN;
        }
        if (component->f <= f_x + xi * beta * gd) {
            evaluate_gradient(g_temp, x_temp, n, component);
            gd_temp = dot_product(g_temp, d, n);
            sigma_gd = sigma * gd;
            if (sigma_gd <= gd_temp) {
                if (-sigma_gd <= fabs(gd_temp)) {
                    component->alpha = beta;
                    free(x_temp);
                    free(g_temp);
                    return LINE_SEARCH_SATISFIED;
                } else {
                    width = decreasing;
                }
            } else {
                width = increasing;
            }
        } else {
            width = decreasing;
        }
        beta *= width;
    }
    free(x_temp);
    free(g_temp);
    return LINE_SEARCH_FAILED;
}

#if 0
/*
 * Name: check_line_search
 * Description:
 * Return:
 *     0 satisfied
 *     2 LINESEARCH_MAXIMUM_STEPWIDTH
 *     3 LINESEARCH_MINIMUM_STEPWIDTH
 */
int
check_line_search(double beta, NonLinearComponent component)
{
    if ( beta < Com->Param->ls_min_stepw ) {
        return LINESEARCH_MINIMUM_STEPWIDTH;
    }
    else if ( beta > Com->Param->ls_max_stepw ) {
        return LINESEARCH_MAXIMUM_STEPWIDTH;
    }
    return 0;
}
#endif
