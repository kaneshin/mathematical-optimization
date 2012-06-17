// vim:set ts=8 sts=4 sw=4 tw=0:
// vim:set foldmethod=marker foldmarker={{{,}}}:
/* ===========================================================================
 *  File: quasi_newton_bfgs
 *  Version: 0.9.0
 *  Last Change: 18-Jun-2012.
 *  Maintainer: Shintaro Kaneko <kaneshin0120@gmail.com>
 *  Description:
=========================================================================== */

#include "include/quasi_newton_bfgs.h"

#include <stdio.h>
#include <stdlib.h>

#include "include/line_search.h"
#include "include/mymath.h"
#include "include/myvector.h"
#include "include/mymatrix.h"

#define DEBUG 0
/*
 * Switch Debug mode
 *      0: OFF
 *      1: ON
 */
#if DEBUG
#   define debug_printf printf
#else
#   define debug_printf 1 ? (void) 0 : printf
#endif

#define _OPTIMIZATION_SEARCH 4
/*
 * Set _OPTIMIZATION_SEARCH
 *      1: armijo
 *      2: wolfe
 *      3: strong wolfe
 *      4: backtracking wolfe
 *      5: backtracking storng wolfe
 */
#if _OPTIMIZATION_SEARCH == 1
#   define _OPTIMIZATION_LINESEARCH \
        line_search_armijo(x, g, d, n, 1., .5, .5, &component)
#elif _OPTIMIZATION_SEARCH == 2
#   define _OPTIMIZATION_LINESEARCH \
        line_search_wolfe(x, g, d, n, 1., .5, .001, .2, &component)
#elif _OPTIMIZATION_SEARCH == 3
#   define _OPTIMIZATION_LINESEARCH \
        line_search_strong_wolfe(x, g, d, n, 1., .5, .001, .2, &component)
#elif _OPTIMIZATION_SEARCH == 4
#   define _OPTIMIZATION_LINESEARCH \
        line_search_backtracking_wolfe(x, g, d, n, 1., .5, .001, .2,\
                parameter->decreasing, parameter->increasing, &component)
#elif _OPTIMIZATION_SEARCH == 5
#   define _OPTIMIZATION_LINESEARCH \
        line_search_backtracking_strong_wolfe(x, g, d, n, 1., .5, .001, .2,\
                parameter->decreasing, parameter->increasing, &component)
#endif

#define _OPTIMIZATION_FORMULA 2
/*
 * Set _OPTIMIZATION_FORMULA
 *      1: B formula
 *      2: H formula
 */
#if _OPTIMIZATION_FORMULA == 1
#   define _OPTIMIZATION_DIRECTION_SEARCH_BFGS \
        direction_search_B_formula(d, b, g, n)
#   define _OPTIMIZATION_UPDATE_BFGS \
        update_bfgs_B_formula(b, s, y, x, n)
#elif _OPTIMIZATION_FORMULA == 2
#   define _OPTIMIZATION_DIRECTION_SEARCH_BFGS \
        direction_search_H_formula(d, b, g, n)
#   define _OPTIMIZATION_UPDATE_BFGS \
        update_bfgs_H_formula(b, s, y, x, n)
#endif

static void quasi_newton_bfgs_parameter_default
(
    QuasiNewtonBFGSParameter *parameter
);

static void initialize_quasi_newton_bfgs
(
    NonLinearComponent *component
);

static int direction_search_B_formula
(
    double *d,
    double **B,
    double *g,
    int n
);

static int update_bfgs_B_formula
(
    double **B,
    const double *s,
    const double *y,
    double *Bs,
    int n
);

static int direction_search_H_formula
(
    double *d,
    double **H,
    const double *g,
    int n
);

static int update_bfgs_H_formula
(
    double **H,
    const double *s,
    const double *y,
    double *Hy,
    int n
);

static void print_iteration_info
(
    int iteration,
    double g_norm,
    NonLinearComponent *component
);

static void print_result_info
(
    int status,
    int iteration,
    NonLinearComponent *component
);

/*
 * quasi-Newton method (BFGS)
 */
int
quasi_newton_bfgs
(
    double *x,
    double **b,
    int n,
    FunctionObject *function_object,
    QuasiNewtonBFGSParameter *parameter
)
{
    /*
     * Return:
     *      QUASI_NEWTON_BFGS_SATISFIED
     *      QUASI_NEWTON_BFGS_OUT_OF_MEMORY
     */
    int i, j, iteration, status;
    double *storage, *storage_x, **storage_b,
           *d, *g, *x_temp, *g_temp, *s, *y,
           g_norm;
    QuasiNewtonBFGSParameter _parameter;
    NonLinearComponent component;

    if (NULL == x) {
        // allocate memory to x as a vector
        if (NULL == (storage_x = (double *)malloc(sizeof(double) * n))) {
            status = QUASI_NEWTON_BFGS_OUT_OF_MEMORY;
            goto result;
        }
        zero_vector(storage_x, n);
        x = storage_x;
    } else {
        storage_x = NULL;
    }
    if (NULL == b) {
        // allocate memory to b as a matrix
        if (NULL == (storage_b = (double **)malloc(sizeof(double *) * n))) {
            status = QUASI_NEWTON_BFGS_OUT_OF_MEMORY;
            goto result;
        }
        if (NULL == (*storage_b = (double *)malloc(sizeof(double) * n * n))) {
            status = QUASI_NEWTON_BFGS_OUT_OF_MEMORY;
            goto result;
        }
        for (i = 1; i < n; i++) {
            storage_b[i] = storage_b[i - 1] + n;
        }
        identity_matrix(storage_b, n);
        b = storage_b;
    } else {
        storage_b = NULL;
    }
    // allocate memory to storage for d, g, xtemp, gtemp, s and y
    if (NULL == (storage = (double *)malloc(sizeof(double) * n * 6))) {
        status = QUASI_NEWTON_BFGS_OUT_OF_MEMORY;
        goto result;
    }
    y = (s = (g_temp = (x_temp = (g = (d = storage) + n) + n) + n) + n) + n;

    if (NULL == parameter) {
        parameter = &_parameter;
        quasi_newton_bfgs_parameter_default(parameter);
    }

    if (NULL == function_object) {
        status = QUASI_NEWTON_BFGS_NON_FUNCTION;
        goto result;
    }
    component.function_object = function_object;

    initialize_quasi_newton_bfgs(&component);

    if (NON_LINEAR_FUNCTION_NAN == evaluate_gradient(g, x, n, &component)) {
        status = QUASI_NEWTON_BFGS_FUNCTION_NAN;
        goto result;
    }
    for (iteration = 1; iteration < parameter->max_iteration; iteration++) {
        // search a direction of descent
        status = _OPTIMIZATION_DIRECTION_SEARCH_BFGS;
        if (status) goto result;
        // compute step width with a line search algorithm
        switch (_OPTIMIZATION_LINESEARCH) {
            case LINE_SEARCH_FUNCTION_NAN:
                status = QUASI_NEWTON_BFGS_FUNCTION_NAN;
                goto result;
            case LINE_SEARCH_FAILED:
                status = QUASI_NEWTON_BFGS_LINE_SEARCH_FAILED;
                goto result;
            default:
                break;
        }
        update_step_vector(x_temp, x, component.alpha, d, n);

        if (NON_LINEAR_FUNCTION_NAN ==
                evaluate_gradient(g_temp, x_temp, n, &component)) {
            status = QUASI_NEWTON_BFGS_FUNCTION_NAN;
            goto result;
        }
        g_norm = norm_infty(g_temp, n);

        print_iteration_info(iteration, g_norm, &component);

        if (g_norm < parameter->tolerance) {
            status = QUASI_NEWTON_BFGS_SATISFIED;
            goto result;
        }

        for (i = 0; i < n; i++) {
            s[i] = x_temp[i] - x[i];
            y[i] = g_temp[i] - g[i];
        }
        switch (_OPTIMIZATION_UPDATE_BFGS) {
            case QUASI_NEWTON_BFGS_FUNCTION_NAN:
                goto result;
            case QUASI_NEWTON_BFGS_SATISFIED:
                printf("* Matrix is NOT updated on BFGS\n");
                break;
            default:
                break;
        }

        copy_vector(x, x_temp, n);
        copy_vector(g, g_temp, n);
    }

result:
    print_result_info(status, iteration, &component);

    if (NULL != storage) {
        free(storage);
        storage = NULL;
    }
    vector_delete(storage_x);
    matrix_delete(storage_b);
    return status;
}

void
quasi_newton_bfgs_parameter_default
(
    QuasiNewtonBFGSParameter *parameter
)
{
    parameter->tolerance = 1.e-9;
    parameter->max_iteration = 1070;
    parameter->step_width = 1.;
    parameter->xi = 0.001;
    parameter->tau = .5;
    parameter->sigma = .2;
    parameter->decreasing = .5;
    parameter->increasing = 2.1;
}

void
initialize_quasi_newton_bfgs
(
    NonLinearComponent *component
)
{
    component->alpha = 0.;
    component->iteration_f  = 0;
    component->iteration_g  = 0;
}

int
direction_search_B_formula
(
    double *d,
    double **B,
    double *g,
    int n
)
{
    /*
     * Search the direction on B formula
     * Return: d: Bd = -g
     */
    int i , status;

    for (i = 0; i < n; i++) g[i] = -g[i];
    status = successive_over_relaxation(B, d, g, n, 1.e-7, 0.5);
    for (i = 0; i < n; i++) g[i] = -g[i];
    return status;
}

int
update_bfgs_B_formula
(
    double **B,
    const double *s,
    const double *y,
    double *Bs,
    int n
)
{
    /*
     * Update the matrix on B formula
     * Return:
     */
    int i, j;
    double temp, sBs, sy;

    for (i = 0; i < n; i++) {
        temp = 0.;
        for (j = 0; j < n; j++) {
            temp += B[i][j] * s[j];
        }
        if (temp != temp) return QUASI_NEWTON_BFGS_FUNCTION_NAN;
        Bs[i] = temp;
    }
    sBs = dot_product(s, Bs, n);
    sy = dot_product(s, y, n);
    if (sy > 0) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                B[i][j] += - Bs[i] * Bs[j] / sBs + y[i] * y[j] / sy;
            }
        }
        return QUASI_NEWTON_BFGS_SATISFIED;
    }
    return QUASI_NEWTON_BFGS_NOT_UPDATE;
}

int
direction_search_H_formula
(
    double *d,
    double **H,
    const double *g,
    int n
)
{
    /*
     * Search the direction on H formula
     * Return: d: d = -Hg
     */
    int i, j;
    double temp;
    for (i = 0; i < n; i++) {
        temp = 0.;
        for (j = 0; j < n; j++) {
            temp -= H[i][j] * g[j];
        }
        if (temp != temp) return QUASI_NEWTON_BFGS_FUNCTION_NAN;
        d[i] = temp;
    }
    return QUASI_NEWTON_BFGS_SATISFIED;
}

int
update_bfgs_H_formula
(
    double **H,
    const double *s,
    const double *y,
    double *Hy,
    int n
)
{
    /*
     * Update the matrix on H formula
     * Return:
     */
    int i, j;
    double temp, yHy, sy;

    for (i = 0; i < n; i++) {
        temp = 0.;
        for (j = 0; j < n; j++) {
            temp += H[i][j] * y[j];
        }
        if (temp != temp) return QUASI_NEWTON_BFGS_FUNCTION_NAN;
        Hy[i] = temp;
    }
    yHy = dot_product(y, Hy, n);
    sy = dot_product(s, y, n);
    if (sy > 0) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                H[i][j] += - (Hy[i] * s[j] + s[i] * Hy[j]) / sy
                                + (1 + yHy / sy) * s[i] * s[j] / sy ;
            }
        }
        return QUASI_NEWTON_BFGS_SATISFIED;
    }
    return QUASI_NEWTON_BFGS_NOT_UPDATE;
}

void
print_iteration_info
(
    int iteration,
    double g_norm,
    NonLinearComponent *component
)
{
    printf("\n == iteration: %7d ================================\n", iteration);
    printf("step width parameter:  \t%13.6e\n", component->alpha);
    printf("-------------------------------------------------------\n");
    printf("function value:        \t%13.6e\n", component->f);
    printf("||gf(x_k+1)||_infinity = %e\n", g_norm);
    printf("-------------------------------------------------------\n");
}

void
print_result_info
(
    int status,
    int iteration,
    NonLinearComponent *component
)
{
    printf("\n == Result ==\n");
    printf("Compute status: %3d\n", status);
    switch (status) {
        case QUASI_NEWTON_BFGS_SATISFIED:
            printf("Satisfied: Quasi-Newton Method (BFGS) is finished\n");
            break;
        case QUASI_NEWTON_BFGS_NON_FUNCTION:
            printf("Failed: Function Object is not defined\n");
            break;
        case QUASI_NEWTON_BFGS_OUT_OF_MEMORY:
            printf("Failed: Out of Memory\n");
            break;
        case QUASI_NEWTON_BFGS_FUNCTION_NAN:
            printf("Failed: Function value is Not a Number\n");
            break;
        case QUASI_NEWTON_BFGS_FAILED:
            printf("Failed: FAILED\n");
            break;
        default:
            break;
    }
    if (status >= QUASI_NEWTON_BFGS_SATISFIED) {
        printf("-------------------------------------------------------\n");
        printf("iterations:          %12d\n", iteration);
        printf("function evaluations:%12d\n", component->iteration_f);
        printf("gradient evaluations:%12d\n", component->iteration_g);
        printf("=======================================================\n");
        printf("function value:      \t%13.6e\n", component->f);
    }
}

