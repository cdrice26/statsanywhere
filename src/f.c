#include "f.h"
#include "special_functions.h"

/**
 * @brief F-distribution probability density function.
 *
 * Evaluates the PDF of the F-distribution with d1 and d2 degrees of
 * freedom at x.
 *
 * @param x Point to evaluate
 * @param d1 Numerator degrees of freedom
 * @param d2 Denominator degrees of freedom
 * @return PDF value
 */
double fpdf(double x, double d1, double d2) {
    return (pow(d1 / d2, d1 / 2) * pow(x, d1 / 2 - 1) * pow(1 + d1 / d2 * x, -(d1 + d2) / 2)) /
    beta(d1 / 2, d2 / 2);
}

/**
 * @brief F-distribution cumulative distribution function (lower tail).
 *
 * Computes the CDF P(X <= x) for X ~ F(d1, d2) using the regularized
 * incomplete beta relationship.
 *
 * @param x Upper limit of integration
 * @param d1 Numerator degrees of freedom
 * @param d2 Denominator degrees of freedom
 * @return CDF value
 */
double fcdf(double x, double d1, double d2) {
    double rbeta_upper_limit = (d1 * x) / (d1 * x + d2);
    return regularized_incomplete_beta(rbeta_upper_limit, d1 / 2.0, d2 / 2.0);
}

int main() {
    printf("(Should be .071589) fpdf(3, 5, 4) = %f\n", fpdf(3, 5, 4));
    printf("(Should be .8453) fcdf(3, 5, 4) = %f\n", fcdf(3, 5, 4));
}
