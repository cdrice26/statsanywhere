#ifndef VECTORS_H
#define VECTORS_H
#include <stdlib.h>
#include <math.h>

double* add(double* v1, double* v2, int n);
double* negate(double* v1, int n);
double* subtract(double* v1, double* v2, int n);
double dot_product(double* v1, double* v2, int n);
double* multiply_scalar(double* v, int n, double k);
double* projection(double* of, double* onto, int n);
double magnitude(double* v, int n);
double* copy_vector(double* v, int n);

#endif
