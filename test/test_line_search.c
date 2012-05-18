#include <CUnit/CUnit.h>
#include <CUnit/Console.h>

#include "../src/include/line_search.h"
#include "../src/include/function_object.h"

#include <stdlib.h>

static int n;
static double *x, *g, *d, **a;

static double function_value(double *x, int n);
static void gradient_vector(double *g, double *x, int n);

void
test_line_search_armijo(void)
{
    int i;
    double result, expect, f_x;
    FunctionObject func;

    expect = .25;
    x[0] = 1.;
    x[1] = 2.;
    gradient_vector(g, x, n);
    d[0] = -2.;
    d[1] = -4.;

    func.function = function_value;
    func.gradient = gradient_vector;

    line_search_armijo(x, g, d, n, .5, &func, &result);

    CU_ASSERT_EQUAL(expect, result);
}

double
function_value(double *x, int n)
{
    return x[0] * x[0] + 2 * x[1] * x[1] * x[1] / 3 + 1;
}

void
gradient_vector(double *g, double *x, int n)
{
    g[0] = 2 * x[0];
    g[1] = 2 * x[1] * x[1];
}

int
main(int argc, char* argv[])
{
    int i;

    n = 2;
    x = (double *)malloc(sizeof(double) * n);
    g = (double *)malloc(sizeof(double) * n);
    d = (double *)malloc(sizeof(double) * n);
    a = (double **)malloc(sizeof(double *) * n);
    *a = (double *)malloc(sizeof(double) * n * n);
    for (i = 1; i < n; i++) a[i] = a[i - 1] + n;

    CU_pSuite testSuite;
    CU_initialize_registry();
    testSuite = CU_add_suite("line_search.c TestSuite", NULL, NULL);

    CU_add_test(testSuite, "line_search_armijo Test", test_line_search_armijo);

    CU_console_run_tests();
    CU_cleanup_registry();

    free(x);
    free(g);
    free(d);
    free(*a);
    free(a);

    return 0;
}

