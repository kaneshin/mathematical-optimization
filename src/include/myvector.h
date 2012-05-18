#ifndef OPTIMIZATION_MYVECTOR_H
#define OPTIMIZATION_MYVECTOR_H

int vector_new_double(double *vector, int n);

void vector_delete(double *vector);

void copy_vector(double *x, const double *y, int n);

double dot_product(const double *x, const double *y, int n);

void zero_vector(double *x, int n);

void update_step_vector
(
    double *x_temp,
    const double *x, double alpha, const double *y, int n
);

double norm_infty(const double *x, int n);

#endif // OPTIMIZATION_MYVECTOR_H

