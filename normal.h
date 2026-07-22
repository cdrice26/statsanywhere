#ifndef NORMAL_H
#define NORMAL_H
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "quadrature_weights.h"

double normalpdf_std(double z);
double normalpdf(double z, double mu, double sigma);
#endif
