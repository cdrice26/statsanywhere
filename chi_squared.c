#include "chi_squared.h"
#include "special_functions.h"

/**
 * @brief Chi-squared probability density function.
 *
 * Evaluates the PDF of the chi-squared distribution with df degrees of
 * freedom at x. Returns 0 for x <= 0.
 *
 * @param x Point to evaluate (must be >= 0)
 * @param df Degrees of freedom
 * @return PDF value
 */
double X2pdf(double x, int df) {
    if (x > 0) {
        return (pow(x, df / 2.0 - 1.0) * exp(-x / 2.0)) / (pow(2.0, df / 2.0) * gamma(df / 2.0));
    } else {
        return 0;
    }
}

/**
 * @brief Chi-squared cumulative distribution function.
 *
 * Computes the lower-tail CDF using the lower incomplete gamma function.
 *
 * @param x Upper limit of integration
 * @param df Degrees of freedom
 * @return CDF value P(X <= x)
 */
double X2cdf(double x, int df) {
    return lower_incomplete_gamma(df / 2.0, x / 2.0) / gamma(df / 2.0);
}

int main() {
    printf("(Should be .073) X2pdf(5, 3) = %f\n", X2pdf(5, 3));
    printf("(Should be .828) X2cdf(5, 3) = %f\n", X2cdf(5, 3));
}
