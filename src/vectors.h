#ifndef VECTORS_H
#define VECTORS_H
#include <stdlib.h>
#include <math.h>

double* add(const double* v1, const double* v2, int n);
double* negate(const double* v1, int n);
double* subtract(const double* v1, const double* v2, int n);
double dot_product(const double* v1, const double* v2, int n);
double* multiply_scalar(const double* v, int n, double k);
double* projection(const double* of, const double* onto, int n);
double magnitude(const double* v, int n);
double* copy_vector(const double* v, int n);

#endif
