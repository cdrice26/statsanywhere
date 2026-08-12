#ifndef POISSON_H
#define POISSON_H
#include "operations.h"
#include <math.h>

double poissonpdf(int k, double lambda);
double poissoncdf(double k, double lambda);

#endif
