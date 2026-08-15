#ifndef POISSON_H
#define POISSON_H
#include "operations.h"
#include <math.h>

double poissonpmf(int k, double lambda);
double poissoncdf(double k, double lambda);

#endif
