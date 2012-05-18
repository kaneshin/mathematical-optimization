#include <CUnit/CUnit.h>
#include <CUnit/Console.h>

#include "../src/include/mymath.h"

#include <stdlib.h>

static int n;
static double *x, *y, **a;

void
test_copy_vector(void)
{
    int i;
    double *result, *expect;

    result = (double *)malloc(n * sizeof(double));
    expect = (double *)malloc(n * sizeof(double));

    for (i = 0; i < n; i++) {
        expect[i] = x[i] = i * 1.;
        result[i] = 0.;
    }

    copy_vector(result, x, n);

    for (i = 0; i < n; i++) {
        CU_ASSERT_EQUAL(expect[i], result[i]);
    }

    free(result);
    free(expect);
}

void
test_dot_product(void)
{
    int i;

    for (i = 0; i < n; i++) {
        x[i] = (i - 3) * 1.;
        y[i] = 1.;
    }
    CU_ASSERT_EQUAL(15., dot_product(x, y, n));

    for (i = 0; i < n; i++) {
        x[i] = i;
        y[i] = 1.;
    }
    CU_ASSERT_EQUAL(45., dot_product(x, y, n));
}

void
test_update_step_vector(void)
{
    int i;
    double alpha, *result, *expect;

    result = (double *)malloc(n * sizeof(double));
    expect = (double *)malloc(n * sizeof(double));

    alpha = 2.5;
    for (i = 0; i < n; i++) {
        x[i] = 1.;
        y[i] = i;
        expect[i] = x[i] + alpha * y[i];
        result[i] = 0.;
    }

    update_step_vector(result, x, alpha, y, n);

    for (i = 0; i < n; i++) {
        CU_ASSERT_EQUAL(expect[i], result[i]);
    }

    free(result);
    free(expect);
}

void
test_norm_infty(void)
{
    int i;

    for (i = 0; i < n; i++) {
        x[i] = 0.;
        y[i] = 0.;
    }
    x[5] = 1.;
    y[1] = -1.;
    CU_ASSERT_EQUAL(1., norm_infty(x, n));
    CU_ASSERT_EQUAL(1., norm_infty(y, n));
}

void
test_identity_matrix(void)
{
    int i, j;
    double **expect;

    expect = (double **)malloc(n * sizeof(double *));
    *expect = (double *)malloc(n * n * sizeof(double));
    for (i = 1; i < n; i++) {
        expect[i] = expect[i - 1] + n;
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            expect[i][j] = 0.;
        }
        expect[i][i] = 1.;
    }
    identity_matrix(a, n);

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            CU_ASSERT_EQUAL(expect[i][j], a[i][j]);
        }
    }

    free(*expect);
    free(expect);
}

void
test_gauss_seidel(void)
{
    int i;
    double *expect;

    expect = (double *)malloc(sizeof(double) * n);
    identity_matrix(a, n);
    for (i = 0; i < n; i++) {
        x[i] = 0.;
        expect[i] = y[i] = i * 1.;
    }
    gauss_seidel(a, x, y, n, 1.e-7);
    for (i = 0; i < n; i++) {
        CU_ASSERT_EQUAL(expect[i], x[i]);
    }
    free(expect);
}

void
test_successive_over_relaxation(void)
{
    int i;
    double *expect;

    expect = (double *)malloc(sizeof(double) * n);
    identity_matrix(a, n);
    for (i = 0; i < n; i++) {
        x[i] = 0.;
        expect[i] = y[i] = i * 1.;
    }
    successive_over_relaxation(a, x, y, n, 1.e-7, 1.);
    for (i = 0; i < n; i++) {
        CU_ASSERT_EQUAL(expect[i], x[i]);
    }
    free(expect);
}

int
main(int argc, char* argv[])
{
    int i;

    n = 10;
    x = (double *)malloc(sizeof(double) * n);
    y = (double *)malloc(sizeof(double) * n);
    a = (double **)malloc(sizeof(double *) * n);
    *a = (double *)malloc(sizeof(double) * n * n);
    for (i = 1; i < n; i++) a[i] = a[i - 1] + n;

    CU_pSuite testSuite;
    CU_initialize_registry();
    testSuite = CU_add_suite("mymath.c TestSuite", NULL, NULL);

    CU_add_test(testSuite, "copy_vector Test", test_copy_vector);
    CU_add_test(testSuite, "dot_product Test", test_dot_product);
    CU_add_test(testSuite, "norm_infty Test", test_norm_infty);
    CU_add_test(testSuite, "identity_matrix Test", test_identity_matrix);
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

