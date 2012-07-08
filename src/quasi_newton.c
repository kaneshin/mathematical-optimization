/*
 * vim:set ts=8 sts=4 sw=4 tw=0:
 *
 * File:        quasi_newton.c
 * Version:     0.2.0
 * Maintainer:  Shintaro Kaneko <kaneshin0120@gmail.com>
 * Last Change: 09-Jul-2012.
 * TODO:
 *  Check elements of parameter
 *  Check each function of function_object
 */

#include "include/quasi_newton.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "include/mymath.h"
#include "include/print_message.h"

static char method_name[64] = "Quasi-Newton";

typedef struct _QuasiNewtonFormula {
    int (*direction_search)(
            double *,
            double **,
            double *,
            int
        );
    int (*update_matrix)(
            double **,
            const double *,
            const double *,
            double *,
            int
        );
} QuasiNewtonFormula;

static void
default_quasi_newton_parameter(
    QuasiNewtonParameter *parameter
);

static void
set_quasi_newton_formula(
    QuasiNewtonFormula *quasi_newton_formula,
    QuasiNewtonParameter *parameter
);

static int
direction_search_bfgs_B_formula(
    double *d,
    double **B,
    double *g,
    int n
);

static int
update_matrix_bfgs_B_formula(
    double **B,
    const double *s,
    const double *y,
    double *Bs,
    int n
);

static int
direction_search_bfgs_H_formula(
    double *d,
    double **H,
    double *g,
    int n
);

static int
update_matrix_bfgs_H_formula(
    double **H,
    const double *s,
    const double *y,
    double *Hy,
    int n
);

int
quasi_newton(
    double *x,
    double **b,
    int n,
    FunctionObject *function_object,
    line_search_t line_search,
    LineSearchParameter *line_search_parameter,
    QuasiNewtonParameter *quasi_newton_parameter
) {
    int i, j, iter, status, storage_b_num, storage_num;
    long int memory_size;
    double g_norm,
           *storage, *storage_x, **storage_b,
           *d, *g, *x_temp, *g_temp, *work, *s, *y;
    NonLinearComponent component;
    QuasiNewtonFormula quasi_newton_formula;
    QuasiNewtonParameter _quasi_newton_parameter;
    EvaluateObject evaluate_object;

    /* memory_size is for doing memcpy */
    memory_size = sizeof(double) * n;
    /* prepare a number of vector for storage_b and storage */
    storage_b_num = n;
    storage_num = 6;
    /*
     * allocate memory to storage
     */
    if (NULL == x) {
        /* allocate memory to storage_x for x as a vector */
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
    if (NULL == b) {
        /* allocate memory to storage_b for b as a matrix
         *      size: storage_b_num * n */
        if (NULL == (storage_b = (double **)malloc(
                        sizeof(double *) * storage_b_num))) {
            status = NON_LINEAR_OUT_OF_MEMORY;
            goto result;
        }
        if (NULL == (*storage_b = (double *)malloc(
                        memory_size * storage_b_num))) {
            status = NON_LINEAR_OUT_OF_MEMORY;
            goto result;
        }
        for (i = 1; i < storage_b_num; ++i) {
            storage_b[i] = storage_b[i - 1] + n;
        }
        /* initialize a matrix as identify */
        for (i = 0; i < storage_b_num; ++i) {
            storage_b[i][i] = 1.;
            for (j = 0; j < i; ++j) {
                storage_b[i][j] = 0.;
            }
            for (j = i + 1; j < n; ++j) {
                storage_b[i][j] = 0.;
            }
        }
        b = storage_b;
    } else {
        storage_b = NULL;
    }
    /* allocate memory to storage for d, g, x_temp, g_temp and
     * work (s and y) */
    if (NULL == (storage = (double *)malloc(memory_size * storage_num))) {
        status = NON_LINEAR_OUT_OF_MEMORY;
        goto result;
    }
    d = storage;
    g = d + n;
    x_temp = g + n;
    g_temp = x_temp + n;
    /* work share memory with s and y */
    s = work = g_temp + n;
    y = s + n;

    /* make sure that f and gf of this problem exist */
    if (NULL == function_object->function
            || NULL == function_object->gradient) {
        status = NON_LINEAR_NO_FUNCTION;
        goto result;
    }

    /* set the component of Non-Linear Programming */
    initialize_non_linear_component(
            method_name, function_object, &evaluate_object, &component);

    /* set the parameter of Quasi-Newton method */
    if (NULL == quasi_newton_parameter) {
        quasi_newton_parameter = &_quasi_newton_parameter;
        default_quasi_newton_parameter(quasi_newton_parameter);
    }

    /* set the formula for solving this problem */
    set_quasi_newton_formula(&quasi_newton_formula, quasi_newton_parameter);

    /* parameter of Line Search */
    if (NULL == line_search_parameter) {
        status = NON_LINEAR_NO_PARAMETER;
        goto result;
    }

    /*
     * start to compute for solving this problem
     */
    if (NON_LINEAR_FUNCTION_OBJECT_NAN
            == evaluate_object.gradient(g, x, n, &component)) {
        status = NON_LINEAR_FUNCTION_NAN;
        goto result;
    }
    for (iter = 1; iter <= quasi_newton_parameter->upper_iter; ++iter) {
        /* search a direction of descent */
        status = quasi_newton_formula.direction_search(d, b, g, n);
        if (status) {
            goto result;
        }
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
        /* compute g_norm */
        g_norm = infinity_norm(g_temp, n);

        print_iteration_info(iter, g_norm, &component);

        if (g_norm < quasi_newton_parameter->tolerance) {
            status = NON_LINEAR_SATISFIED;
            goto result;
        }

        /* compute s = x_temp - x and y = g_temp - g */
        for (i = 0; i < n; ++i) {
            s[i] = x_temp[i] - x[i];
            y[i] = g_temp[i] - g[i];
        }
        /* update matrix */
        status = quasi_newton_formula.update_matrix(b, s, y, x, n);
        switch (status) {
            case NON_LINEAR_FUNCTION_NAN:
                goto result;
            case NON_LINEAR_NOT_UPDATE:
                printf("* Matrix is NOT updated\n");
            case NON_LINEAR_SATISFIED:
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
default_quasi_newton_parameter(
    QuasiNewtonParameter *parameter
) {
    parameter->formula = 'h';
    parameter->tolerance = 1.e-8;
    parameter->upper_iter = 5000;
}

static void
set_quasi_newton_formula(
    QuasiNewtonFormula *quasi_newton_formula,
    QuasiNewtonParameter *parameter
) {
    switch (parameter->formula) {
        case 'b': case 'B':
            quasi_newton_formula->direction_search = direction_search_bfgs_B_formula;
            quasi_newton_formula->update_matrix = update_matrix_bfgs_B_formula;
            break;
        case 'h': case 'H':
            quasi_newton_formula->direction_search = direction_search_bfgs_H_formula;
            quasi_newton_formula->update_matrix = update_matrix_bfgs_H_formula;
            break;
        default:
            quasi_newton_formula->direction_search = direction_search_bfgs_H_formula;
            quasi_newton_formula->update_matrix = update_matrix_bfgs_H_formula;
            break;
    }
}

static int
direction_search_bfgs_B_formula(
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
update_matrix_bfgs_B_formula(
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
        if (Bsi != Bsi) return NON_LINEAR_FUNCTION_NAN;
        Bs[i] = Bsi;
    }
    for (i = 1, sBs = s[0] * Bs[0], sy = s[0] * y[0]; i < n; ++i) {
        sBs += s[i] * Bs[i];
        sy += s[i] * y[i];
    }
    if (sBs != sBs || sy != sy) return NON_LINEAR_FUNCTION_NAN;
    if (sy > 0) {
        for (i = 0; i < n; ++i) {
            for (j = 0; j < n; j++) {
                B[i][j] += - Bs[i] * Bs[j] / sBs + y[i] * y[j] / sy;
            }
        }
        return NON_LINEAR_SATISFIED;
    }
    return NON_LINEAR_NOT_UPDATE;
}

static int
direction_search_bfgs_H_formula(
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
        if (Hg != Hg) return NON_LINEAR_FUNCTION_NAN;
        d[i] = Hg;
    }
    return NON_LINEAR_SATISFIED;
}

static int
update_matrix_bfgs_H_formula(
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
        if (Hyi != Hyi) return NON_LINEAR_FUNCTION_NAN;
        Hy[i] = Hyi;
    }
    for (i = 1, yHy = y[0] * Hy[0], sy = s[0] * y[0]; i < n; ++i) {
        yHy += y[i] * Hy[i];
        sy += s[i] * y[i];
    }
    if (yHy != yHy || sy != sy) return NON_LINEAR_FUNCTION_NAN;
    if (sy > 0) {
        for (i = 0; i < n; ++i) {
            for (j = 0; j < n; j++) {
                H[i][j] += - (Hy[i] * s[j] + s[i] * Hy[j]) / sy
                                + (1 + yHy / sy) * s[i] * s[j] / sy ;
            }
        }
        return NON_LINEAR_SATISFIED;
    }
    return NON_LINEAR_NOT_UPDATE;
}

