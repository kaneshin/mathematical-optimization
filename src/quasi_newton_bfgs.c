/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        quasi_newton_bfgs.c
 * Version:     0.1.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 02-Jul-2012.
 * TODO:
 *  Check elements of parameter
 *  Check each function of function_object
 */

#include "include/quasi_newton_bfgs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "include/mymath.h"

typedef struct _BFGSFormula {
    int (*direction_search)(
            double *,
            double **,
            double *,
            int
        );
    int (*update_bfgs)(
            double **,
            const double *,
            const double *,
            double *,
            int
        );
} BFGSFormula;

static void
default_quasi_newton_bfgs_parameter(
    QuasiNewtonBFGSParameter *parameter
);

static void
set_bfgs_formula(
    BFGSFormula *bfgs_formula,
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
    unsigned int n,
    FunctionObject *function_object,
    line_search_t line_search,
    LineSearchParameter *line_search_parameter,
    QuasiNewtonBFGSParameter *quasi_newton_bfgs_parameter
) {
    int status;
    unsigned int i, j, iter;
    unsigned long int memory_size;
    double *storage, *storage_x, **storage_b,
           *d, *g, *x_temp, *g_temp, *s, *y, g_norm, *work;
    NonLinearComponent component;
    BFGSFormula bfgs_formula;
    QuasiNewtonBFGSParameter _quasi_newton_bfgs_parameter;
    EvaluateObject evaluate_object;

    memory_size = sizeof(double) * n;
    /* allocate memory to storage for x_temp and g_temp */
    if (NULL == x) {
        /* allocate memory to storage_x */
        if (NULL == (storage_x = (double *)malloc(memory_size))) {
            status = QUASI_NEWTON_BFGS_OUT_OF_MEMORY;
            goto result;
        }
        for (i = 0; i < n; ++i) {
            storage_x[i] = 0.;
        }
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
        if (NULL == (*storage_b = (double *)malloc(memory_size * n))) {
            status = QUASI_NEWTON_BFGS_OUT_OF_MEMORY;
            goto result;
        }
        for (i = 1; i < n; ++i) {
            storage_b[i] = storage_b[i - 1] + n;
        }
        for (i = 0; i < n; ++i) {
            for (j = 0; j < n; ++j) {
                storage_b[i][j] = 0.;
            }
        }
        b = storage_b;
    } else {
        storage_b = NULL;
    }
    /* allocate memory to storage for d, g, x_temp, g_temp, s, y and work */
    if (NULL == (storage = (double *)malloc(memory_size * 8))) {
        status = QUASI_NEWTON_BFGS_OUT_OF_MEMORY;
        goto result;
    }
    d = storage;
    g = d + n;
    x_temp = g + n;
    g_temp = x_temp + n;
    s = g_temp + n;
    y = s + n;
    work = y + n;

    /* make sure that f and gf of this problem exist */
    if (NULL == function_object->function || NULL == function_object->gradient) {
        status = QUASI_NEWTON_BFGS_NO_FUNCTION;
        goto result;
    }

    /* set the component of Non-Linear Programming */
    initialize_non_linear_component(function_object, &evaluate_object, &component);

    /* TODO: Set defaule parameter in default_quasi_newton_bfgs_parameter
     * set the parameter of Quasi-Newton method on BFGS */
    if (NULL == quasi_newton_bfgs_parameter) {
        quasi_newton_bfgs_parameter = &_quasi_newton_bfgs_parameter;
        default_quasi_newton_bfgs_parameter(quasi_newton_bfgs_parameter);
    }

    /* set the formula for solving this problem */
    set_bfgs_formula(&bfgs_formula, quasi_newton_bfgs_parameter);

    /* parameter of Line Search */
    if (NULL == line_search_parameter) {
        status = QUASI_NEWTON_BFGS_NO_PARAMETER;
        goto result;
    }

    /* start to compute for solving this problem */
    if (NON_LINEAR_FUNCTION_NAN == evaluate_object.gradient(g, x, n, &component)) {
        status = QUASI_NEWTON_BFGS_FUNCTION_NAN;
        goto result;
    }
    for (iter = 1; iter <= quasi_newton_bfgs_parameter->upper_iter; ++iter) {
        /* search a direction of descent */
        if (status = bfgs_formula.direction_search(d, b, g, n)) {
            goto result;
        }
        /* compute step width with a line search algorithm */
        switch (line_search(work, x, g, d, n,
                    &evaluate_object, line_search_parameter, &component)) {
            case LINE_SEARCH_FUNCTION_NAN:
                status = QUASI_NEWTON_BFGS_FUNCTION_NAN;
                goto result;
            case LINE_SEARCH_FAILED:
                status = QUASI_NEWTON_BFGS_LINE_SEARCH_FAILED;
                goto result;
            default:
                break;
        }
        /* update x_temp = x + alpha * d */
        for (i = 0; i < n; ++i) {
            x_temp[i] = x[i] + component.alpha * d[i];
        }
        /* update g_temp = gradient(x_temp) */
        if (NON_LINEAR_FUNCTION_NAN ==
                evaluate_object.gradient(g_temp, x_temp, n, &component)) {
            status = QUASI_NEWTON_BFGS_FUNCTION_NAN;
            goto result;
        }
        /*TODO: compute g_norm */
        g_norm = infinity_norm(g_temp, n);

        print_iteration_info(iter, g_norm, &component);

        if (g_norm < quasi_newton_bfgs_parameter->tolerance) {
            status = QUASI_NEWTON_BFGS_SATISFIED;
            goto result;
        }

        /* compute s = x_temp - x and y = g_temp - g */
        for (i = 0; i < n; ++i) {
            s[i] = x_temp[i] - x[i];
            y[i] = g_temp[i] - g[i];
        }
        /* update matrix of bfgs */
        switch (status = bfgs_formula.update_bfgs(b, s, y, x, n)) {
            case QUASI_NEWTON_BFGS_FUNCTION_NAN:
                goto result;
            case QUASI_NEWTON_BFGS_NOT_UPDATE:
                printf("* Matrix is NOT updated on BFGS\n");
            case QUASI_NEWTON_BFGS_SATISFIED:
            default:
                break;
        }

        memcpy(x, x_temp, memory_size);
        memcpy(g, g_temp, memory_size);
    }

result:
    print_result_info(status, iter, &component);

    /* release memory of storage_x, storage_b and storage_x */
    if (NULL != storage_x) {
        free(storage_x);
        storage_x = NULL;
    }
    if (NULL != storage_b) {
        if (NULL != *storage_b) {
            free(*storage_b);
            *storage_b = NULL;
        }
        free(storage_b);
        storage_b = NULL;
    }
    if (NULL != storage) {
        free(storage);
        storage = NULL;
    }

    return status;
}

static void
default_quasi_newton_bfgs_parameter(
    QuasiNewtonBFGSParameter *parameter
) {
    parameter->formula = 'h';
    parameter->tolerance = 1.e-8;
    parameter->upper_iter = 5000;
}

static void
set_bfgs_formula(
    BFGSFormula *bfgs_formula,
    QuasiNewtonBFGSParameter *parameter
) {
    switch (parameter->formula) {
        case 'b': case 'B':
            bfgs_formula->direction_search = direction_search_B_formula;
            bfgs_formula->update_bfgs = update_bfgs_B_formula;
            break;
        case 'h': case 'H':
            bfgs_formula->direction_search = direction_search_H_formula;
            bfgs_formula->update_bfgs = update_bfgs_H_formula;
            break;
        default:
            bfgs_formula->direction_search = direction_search_H_formula;
            bfgs_formula->update_bfgs = update_bfgs_H_formula;
            break;
    }
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
        case QUASI_NEWTON_BFGS_NO_PARAMETER:
            printf("Failed: Parameter of line search is not defined\n");
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

