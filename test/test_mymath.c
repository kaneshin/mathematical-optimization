#include <CUnit/CUnit.h>
#include <CUnit/Console.h>

#include "../src/include/mymath.h"

#include <stdlib.h>

static int n;
static double *x, *y, **a;

void
test_dot_product(void) {
    /* double
     * dot_product(
     *     const double *x,
     *     const double *y,
     *     unsigned int n
     * ); */
    unsigned int i;

    for (i = 0; i < n; ++i) {
        x[i] = (i - 3) * 1.;
        y[i] = 1.;
    }
    printf("\nkanesihn: %f\n", dot_product(x, y, n));
    CU_ASSERT_EQUAL(15., dot_product(x, y, n));

    for (i = 0; i < n; ++i) {
        y[i] = i;
        x[i] = 1.;
    }
    CU_ASSERT_EQUAL(45., dot_product(x, y, n));
}

void
test_update_step_vector(void) {
    /* void
     * update_step_vector(
     *     double *x_temp,
     *     const double *x,
     *     double alpha,
     *     const double *y,
     *     unsigned int n
     * ); */
    unsigned int i;
    double alpha, *result, *expect;

    result = (double *)malloc(n * sizeof(double));
    expect = (double *)malloc(n * sizeof(double));

    alpha = 2.5;
    for (i = 0; i < n; ++i) {
        x[i] = 1.;
        y[i] = i;
        expect[i] = x[i] + alpha * y[i];
        result[i] = 0.;
    }
    update_step_vector(result, x, alpha, y, n);
    for (i = 0; i < n; ++i) {
        CU_ASSERT_EQUAL(expect[i], result[i]);
    }

    free(result);
    free(expect);
}

void
test_manhattan_norm(void) {
    /* double
     * manhattan_norm(
     *     const double *x,
     *     unsigned int n
     * ); */
    int i;

    for (i = 0; i < n; ++i) {
        x[i] = 0.;
        y[i] = 0.;
    }
    x[5] = 1.;
    y[1] = -2.; y[2] = -4.;
    CU_ASSERT_EQUAL(1., manhattan_norm(x, n));
    CU_ASSERT_EQUAL(6., manhattan_norm(y, n));
}

void
test_euclidean_norm(void) {
    /* double
     * euclidean_norm(
     *     const double *x,
     *     unsigned int n
     * ); */
    int i;

    for (i = 0; i < n; ++i) {
        x[i] = 0.;
        y[i] = 0.;
    }
    x[5] = 1.;
    y[1] = -2.;
    CU_ASSERT_EQUAL(1., euclidean_norm(x, n));
    CU_ASSERT_EQUAL(2., euclidean_norm(y, n));
}

void
test_infinity_norm(void) {
    /* double
     * infinity_norm(
     *     const double *x,
     *     unsigned int n
     * ); */
    int i;

    for (i = 0; i < n; ++i) {
        x[i] = 0.;
        y[i] = 0.;
    }
    x[5] = 1.;
    y[1] = -2.;
    CU_ASSERT_EQUAL(1., infinity_norm(x, n));
    CU_ASSERT_EQUAL(2., infinity_norm(y, n));
}

void
test_gauss_seidel(void) {
    /* int
     * gauss_seidel(
     * double **a,
     * double *x,
     * const double *b,
     * unsigned int n,
     * double epsilon
     * ); */
    unsigned int i, j;
    double *expect;

    expect = (double *)malloc(sizeof(double) * n);
    for (i = 0; i < n; ++i) {
        a[i][i] = 1.;
        for (j = 0; j < i; ++j) {
            a[i][j] = 0.;
        }
        for (j = i + 1; j < n; ++j) {
            a[i][j] = 0.;
        }
    }
    for (i = 0; i < n; ++i) {
        x[i] = 0.;
        expect[i] = y[i] = i * 1.;
    }
    gauss_seidel(a, x, y, n, 1.e-7);
    for (i = 0; i < n; ++i) {
        CU_ASSERT_EQUAL(expect[i], x[i]);
    }
    free(expect);
}

void
test_successive_over_relaxation(void) {
    /* int
     * successive_over_relaxation(
     *     double **a,
     *     double *x,
     *     const double *b,
     *     unsigned int n,
     *     double epsilon,
     *     double omega
     * ); */
    unsigned int i, j;
    double *expect;

    expect = (double *)malloc(sizeof(double) * n);
    for (i = 0; i < n; ++i) {
        a[i][i] = 1.;
        for (j = 0; j < i; ++j) {
            a[i][j] = 0.;
        }
        for (j = i + 1; j < n; ++j) {
            a[i][j] = 0.;
        }
    }
    for (i = 0; i < n; ++i) {
        x[i] = 0.;
        expect[i] = y[i] = i * 1.;
    }
    successive_over_relaxation(a, x, y, n, 1.e-7, 1.);
    for (i = 0; i < n; ++i) {
        CU_ASSERT_EQUAL(expect[i], x[i]);
    }
    free(expect);
}

int
main(int argc, char* argv[]) {
    int i;

    n = 10;
    x = (double *)malloc(sizeof(double) * n);
    y = (double *)malloc(sizeof(double) * n);
    a = (double **)malloc(sizeof(double *) * n);
    *a = (double *)malloc(sizeof(double) * n * n);
    for (i = 1; i < n; ++i) a[i] = a[i - 1] + n;

    CU_pSuite testSuite;
    CU_initialize_registry();
    testSuite = CU_add_suite("mymath.c TestSuite", NULL, NULL);

    CU_add_test(testSuite, "dot_product Test", test_dot_product);
    CU_add_test(testSuite, "update_step_vector Test", test_update_step_vector);
    CU_add_test(testSuite, "manhattan_norm Test", test_manhattan_norm);
    CU_add_test(testSuite, "euclidean_norm Test", test_euclidean_norm);
    CU_add_test(testSuite, "infinity_norm Test", test_infinity_norm);
    CU_add_test(testSuite, "gauss_seidel Test", test_gauss_seidel);
    CU_add_test(testSuite, "successive_over_relaxation Test", test_successive_over_relaxation);

    CU_console_run_tests();
    CU_cleanup_registry();

    free(x);
    free(y);
    free(*a);
    free(a);

    return 0;
}

