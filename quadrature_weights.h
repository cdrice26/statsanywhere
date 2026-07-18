#ifndef QUADRATURE_WEIGHTS_H
#define QUADRATURE_WEIGHTS_H
#include <math.h>
#include <stdlib.h>
#include "matrices.h"
double quadrature_estimate(int n, double* a, double* b, double beta_0, double (*f)(double, void*), void* ctx);
#endif
