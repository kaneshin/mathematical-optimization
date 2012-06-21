/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        quasi_newton_bfgs.c
 * Version:     0.1.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 22-Jun-2012.
 * TODO:
 *  Check elements of parameter
 *  Check each function of function_object
 */

#include "include/quasi_newton_bfgs.h"

#include <stdio.h>
#include <stdlib.h>

#include "include/backtracking_wolfe.h"
#include "include/mymath.h"
#include "include/myvector.h"
#include "include/mymatrix.h"

static void
default_quasi_newton_bfgs_parameter(
    QuasiNewtonBFGSParameter *parameter
);

static int
direction_search_B_formula(
    double *d,
    double **B,
    double *g,
    int n
);

static int
update_bfgs_B_formula(
    double **B,
    const double *s,
    const double *y,
    double *Bs,
    int n
);

static int
direction_search_H_formula(
    double *d,
    double **H,
    double *g,
    int n
);

static int
update_bfgs_H_formula(
    double **H,
    const double *s,
    const double *y,
    double *Hy,
    int n
);

static void
print_iteration_info(
    int iteration,
    double g_norm,
    NonLinearComponent *component
);

static void
print_result_info(
    int status,
    int iteration,
    NonLinearComponent *component
);

int
quasi_newton_bfgs(
    double *x,
    double **b,
    int n,
    FunctionObject *function_object,
    QuasiNewtonBFGSParameter *quasi_newton_bfgs_parameter,
    char formula,
    LineSearchParameter *line_search_parameter,
    int (*line_search)(double *, double *, double *, int,
        LineSearchParameter *, NonLinearComponent *
    )
) {
    int i, k, status;
    double *storage, *storage_x, **storage_b,
           *d, *g, *x_temp, *g_temp, *s, *y, g_norm;
    NonLinearComponent component;
    QuasiNewtonBFGSParameter _quasi_newton_bfgs_parameter;
    BFGSFormula bfgs_formula;
    LineSearchParameter _line_search_parameter;

    if (NULL == x) {
        /* allocate memory to storage_x */
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
        /* allocate memory to storage_b */
        if (NULL == (storage_b = (double **)malloc(sizeof(double *) * n))) {
            status = QUASI_NEWTON_BFGS_OUT_OF_MEMORY;
            goto result;
        }
        if (NULL == (*storage_b = (double *)malloc(sizeof(double) * n * n))) {
            status = QUASI_NEWTON_BFGS_OUT_OF_MEMORY;
            goto result;
        }
        for (i = 1; i < n; ++i) {
            storage_b[i] = storage_b[i - 1] + n;
        }
        identity_matrix(storage_b, n);
        b = storage_b;
    } else {
        storage_b = NULL;
    }
    /* allocate memory to storage for d, g, x_temp, g_temp, s and y */
    if (NULL == (storage = (double *)malloc(sizeof(double) * n * 6))) {
        status = QUASI_NEWTON_BFGS_OUT_OF_MEMORY;
        goto result;
    }
    d = storage;
    g = d + n;
    x_temp = g + n;
    g_temp = x_temp + n;
    s = g_temp + n;
    y = s + n;

    /* make sure that f and gf of this problem exist */
    if (NULL == function_object) {
        status = QUASI_NEWTON_BFGS_NO_FUNCTION;
        goto result;
    }
    /* set the component of Non-Linear Programming */
    initialize_non_linear_component(function_object, &component);
    /* set the parameter of Quasi-Newton method on BFGS */
    if (NULL == quasi_newton_bfgs_parameter) {
        quasi_newton_bfgs_parameter = &_quasi_newton_bfgs_parameter;
        default_quasi_newton_bfgs_parameter(quasi_newton_bfgs_parameter);
    }
    /* set the formula for solving this problem */
    switch (formula) {
        case 'b': case 'B':
            bfgs_formula.direction_search = direction_search_B_formula;
            bfgs_formula.update_bfgs = update_bfgs_B_formula;
            break;
        case 'h': case 'H':
            bfgs_formula.direction_search = direction_search_H_formula;
            bfgs_formula.update_bfgs = update_bfgs_H_formula;
            break;
        default:
            bfgs_formula.direction_search = direction_search_H_formula;
            bfgs_formula.update_bfgs = update_bfgs_H_formula;
            break;
    }
    /* parameter of Line Search */
    if (NULL == line_search_parameter) {
        line_search_parameter = &_line_search_parameter;
        default_line_search_parameter(line_search_parameter);
    }

    /* start to compute for solving this problem */
    if (NON_LINEAR_FUNCTION_NAN == evaluate_gradient(g, x, n, &component)) {
        status = QUASI_NEWTON_BFGS_FUNCTION_NAN;
        goto result;
    }
    for (k = 1; k < quasi_newton_bfgs_parameter->upper_iteration; ++k) {
        /* search a direction of descent */
        status = bfgs_formula.direction_search(d, b, g, n);
        if (status) goto result;
        /* compute step width with a line search algorithm */
        switch (line_search(x, g, d, n, line_search_parameter, &component)) {
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

        print_iteration_info(k, g_norm, &component);

        if (g_norm < quasi_newton_bfgs_parameter->tolerance) {
            status = QUASI_NEWTON_BFGS_SATISFIED;
            goto result;
        }

        for (i = 0; i < n; ++i) {
            s[i] = x_temp[i] - x[i];
            y[i] = g_temp[i] - g[i];
        }
        status = bfgs_formula.update_bfgs(b, s, y, x, n);
        switch (status) {
            case QUASI_NEWTON_BFGS_FUNCTION_NAN:
                goto result;
            case QUASI_NEWTON_BFGS_NOT_UPDATE:
                printf("* Matrix is NOT updated on BFGS\n");
            case QUASI_NEWTON_BFGS_SATISFIED:
            default:
                break;
        }

        for (i = 0; i < n; ++i) {
            x[i] = x_temp[i];
            g[i] = g_temp[i];
        }
    }

result:
    print_result_info(status, k, &component);

    if (NULL != storage) {
        free(storage);
        storage = NULL;
    }
    if (NULL != storage_x) {
        free(storage_x);
        storage_x = NULL;
    }
    if (NULL != storage_b) {
        if (NULL != *storage_b) {
            free(*storage_b);
        }
        free(storage_b);
        storage_b = NULL;
    }

    return status;
}

static void
default_quasi_newton_bfgs_parameter(
    QuasiNewtonBFGSParameter *parameter
) {
    parameter->tolerance        = 1.e-7;
    parameter->upper_iteration  = 3000;
}

static int
direction_search_B_formula(
    double *d,
    double **B,
    double *g,
    int n
) {
    int i , status;

    for (i = 0; i < n; ++i) g[i] = -g[i];
    status = successive_over_relaxation(B, d, g, n, 1.e-7, 0.5);
    for (i = 0; i < n; ++i) g[i] = -g[i];
    return status;
}

static int
update_bfgs_B_formula(
    double **B,
    const double *s,
    const double *y,
    double *Bs,
    int n
) {
    int i, j;
    double Bsi, sBs, sy;

    for (i = 0; i < n; ++i) {
        for (j = 0, Bsi = 0.; j < n; j++) {
            Bsi += B[i][j] * s[j];
        }
        if (Bsi != Bsi) return QUASI_NEWTON_BFGS_FUNCTION_NAN;
        Bs[i] = Bsi;
    }
    for (i = 1, sBs = s[0] * Bs[0], sy = s[0] * y[0]; i < n; ++i) {
        sBs += s[i] * Bs[i];
        sy += s[i] * y[i];
    }
    if (sBs != sBs || sy != sy) return QUASI_NEWTON_BFGS_FUNCTION_NAN;
    if (sy > 0) {
        for (i = 0; i < n; ++i) {
            for (j = 0; j < n; j++) {
                B[i][j] += - Bs[i] * Bs[j] / sBs + y[i] * y[j] / sy;
            }
        }
        return QUASI_NEWTON_BFGS_SATISFIED;
    }
    return QUASI_NEWTON_BFGS_NOT_UPDATE;
}

static int
direction_search_H_formula(
    double *d,
    double **H,
    double *g,
    int n
) {
    int i, j;
    double Hg;
    for (i = 0; i < n; ++i) {
        for (j = 0, Hg = 0.; j < n; j++) {
            Hg -= H[i][j] * g[j];
        }
        if (Hg != Hg) return QUASI_NEWTON_BFGS_FUNCTION_NAN;
        d[i] = Hg;
    }
    return QUASI_NEWTON_BFGS_SATISFIED;
}

static int
update_bfgs_H_formula(
    double **H,
    const double *s,
    const double *y,
    double *Hy,
    int n
) {
    int i, j;
    double Hyi, yHy, sy;

    for (i = 0; i < n; ++i) {
        for (j = 0, Hyi = 0.; j < n; j++) {
            Hyi += H[i][j] * y[j];
        }
        if (Hyi != Hyi) return QUASI_NEWTON_BFGS_FUNCTION_NAN;
        Hy[i] = Hyi;
    }
    for (i = 1, yHy = y[0] * Hy[0], sy = s[0] * y[0]; i < n; ++i) {
        yHy += y[i] * Hy[i];
        sy += s[i] * y[i];
    }
    if (yHy != yHy || sy != sy) return QUASI_NEWTON_BFGS_FUNCTION_NAN;
    if (sy > 0) {
        for (i = 0; i < n; ++i) {
            for (j = 0; j < n; j++) {
                H[i][j] += - (Hy[i] * s[j] + s[i] * Hy[j]) / sy
                                + (1 + yHy / sy) * s[i] * s[j] / sy ;
            }
        }
        return QUASI_NEWTON_BFGS_SATISFIED;
    }
    return QUASI_NEWTON_BFGS_NOT_UPDATE;
}

static void
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

static void
print_result_info(
    int status,
    int iteration,
    NonLinearComponent *component
) {
    printf("\n\n\nCompute status: %3d\n", status);
    switch (status) {
        case QUASI_NEWTON_BFGS_SATISFIED:
            printf("Satisfied: Quasi-Newton Method (BFGS) is finished\n");
            break;
        case QUASI_NEWTON_BFGS_FUNCTION_NAN:
            printf("Failed: Function value is Not a Number\n");
            break;
        case QUASI_NEWTON_BFGS_OUT_OF_MEMORY:
            printf("Failed: Out of Memory\n");
            break;
        case QUASI_NEWTON_BFGS_NO_FUNCTION:
            printf("Failed: Function Object is not defined\n");
            break;
        case QUASI_NEWTON_BFGS_FAILED:
            printf("Failed: FAILED\n");
            break;
        case QUASI_NEWTON_BFGS_LINE_SEARCH_FAILED:
            printf("Failed: Line Search is failed\n");
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

