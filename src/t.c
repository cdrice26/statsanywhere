#include "t.h"
#include "special_functions.h"

/**
 * @brief Student's t-distribution probability density function.
 *
 * Evaluates the PDF of the t-distribution with df degrees of freedom at t.
 *
 * @param t Point to evaluate
 * @param df Degrees of freedom
 * @return PDF value
 */
double tpdf(double t, int df) {
    return gamma((df + 1.0) / 2.0) / (sqrt(M_PI * df) * gamma(df / 2.0)) * pow(1 + pow(t, 2.0) / df, -(df + 1.0) / 2.0);
}

/**
 * @brief Student's t-distribution cumulative distribution function.
 *
 * Computes the lower-tail CDF P(T <= t) using a relationship with the
 * regularized incomplete beta function.
 *
 * @warning Not quite as precise as some other functions due to multiple quadratures being used
 *
 * @param t Point to evaluate
 * @param df Degrees of freedom
 * @return CDF value
 */
double tcdf(double t, int df) {
    double x = df / (pow(t, 2) + df);
    return 1.0 - 0.5 * regularized_incomplete_beta(x, df / 2.0, 0.5);
}

int main() {
    printf("(Should be .0173) tpdf(3, 5) = %f\n", tpdf(3, 5));
    printf("(Should be .985) tcdf(3, 5) = %f\n", tcdf(3, 5));
    printf("(Should be .681) tcdf(0.5, 5) = %f\n", tcdf(0.5, 5));
}
