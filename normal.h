#ifndef NORMAL_H
#define NORMAL_H
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "quadrature_weights.h"
#include "special_functions.h"

double normalpdf_std(double z);
double normalpdf(double z, double mu, double sigma);

double normalcdf_std(double lower, double upper);
double normalcdf(double lower, double upper, double mu, double sigma);
#endif
