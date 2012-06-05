#include "driver.h"

#include <stdlib.h>
#include <math.h>

#include "src/include/quasi_newton_bfgs.h"
#include "src/include/non_linear_component.h"
#include "src/include/myvector.h"

static double function(const double *x, int n);
static void gradient(double *g, const double *x, int n);

int
main(int argc, char* argv[])
{
    int i, n;
    double *x;
    FunctionObject Function;

    n = 100;
    x = (double *)malloc(sizeof(double) * n);

    for (i = 0; i < n; i++) x[i] = 1.;

    Function.function = function;
    Function.gradient = gradient;

    quasi_newton_bfgs(x, NULL, n, &Function, NULL);

    vector_delete(x);

    return 0;
}

double
function(const double *x, int n)
{
    int i;
    double f = 0.;
    for (i = 0; i < n; i++) {
        f += exp(x[i]) - x[i] * sqrt(i + 1.);
    }
    return f;
}

void
gradient(double *g, const double *x, int n)
{
    int i;
    for (i = 0; i < n; i++) {
        g [i] = exp(x [i]) - sqrt (i + 1.);
    }
}

