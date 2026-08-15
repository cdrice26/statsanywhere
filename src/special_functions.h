#ifndef SPECIAL_FUNCTIONS_H
#define SPECIAL_FUNCTIONS_H
#include <math.h>
#include "quadrature_weights.h"

double gamma(double z);
double lower_incomplete_gamma(double s, double x);
double beta(double z1, double z2);
double incomplete_beta(double x, double a, double b);
double regularized_incomplete_beta(double x, double a, double b);
double erf(double z);

#endif
