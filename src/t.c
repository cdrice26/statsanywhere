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
    /* Use identity I_x(a,b) = 1 - I_{1-x}(b,a) to improve numerical stability.
       For t >= 0: CDF = 1 - 0.5 * I_x(a,0.5) = 0.5 + 0.5 * I_{1-x}(0.5,a)
       For t < 0: CDF = 0.5 - 0.5 * I_{1-x}(0.5,a)
    */
    double x = df / (pow(t, 2) + df);
    double a = df / 2.0;
    double b = 0.5;
    double ix_comp = regularized_incomplete_beta(1.0 - x, b, a);
    if (t >= 0) {
        return 0.5 + 0.5 * ix_comp;
    } else {
        return 0.5 - 0.5 * ix_comp;
    }
}

