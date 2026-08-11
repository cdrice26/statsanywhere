#ifndef BINOMIAL_H
#define BINOMIAL_H
#include "operations.h"
#include <math.h>

double binompdf(int n, int k, double p);
double binomcdf(int n, double k, double p);

#endif
