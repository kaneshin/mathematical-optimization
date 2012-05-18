#include "driver.h"

#include <stdlib.h>
#include <math.h>

#include "src/include/quasi_newton_bfgs.h"
#include "src/include/non_linear_component.h"
#include "src/include/myvector.h"
#include "src/include/mymatrix.h"

/*
 * Problem:
 *  minimize f(x) = (x1 - x2^2)^2 / 2 + (x2 - 2)^2 / 2
 *
 *  gradient f(x) = [ gf(x)1, gf(x)2 ]^T
 *      gf(x)1 = x1 - x2^2
 *      gf(x)2 = -2x2(x1 - x2^2) + x2 - 2
 *
 *  x* = [ 4, 2 ]^T
 *
 * NOTE:
 *  x := (x1, x2)
 *  xN : elements of x
 *  x* : optimal solution
 */

static double function(const double *x, int n);
static void gradient(double *g, const double *x, int n);

int
main(int argc, char* argv[])
{
    int i, n;
    double *x, **b;

    FunctionObject Function;

    n = 2;
    Function.function = function;
    Function.gradient = gradient;

    x = (double *)malloc(sizeof(double) * n);
    b = (double **)malloc(sizeof(double *) * n);
    *b = (double *)malloc(sizeof(double) * n * n);
    for (i = 1; i < n; i++) b[i] = b[i - 1] + n;

    for (i = 0; i < n; i++) x[i] = 0.;
    b[0][0] = 1.;   b[0][1] = -2.;
    b[1][0] = -2.;  b[1][1] = 6.;

    quasi_newton_bfgs(x, NULL, n, &Function, NULL);

    vector_delete(x);
    matrix_delete(b);

    return 0;
}

double
function(const double *x, int n)
{
    return (pow((x[0] - pow(x[1], 2)), 2) + pow(x[1] - 2, 2)) / 2.;
}

void
gradient(double *g, const double *x, int n)
{
    g[0] = x[0] - pow(x[1], 2);
    g[1] = -2. * x[1] * (x[0] - pow(x[1], 2)) + x[1] - 2.;
}

