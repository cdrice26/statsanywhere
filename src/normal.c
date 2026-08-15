#include "normal.h"

/**
 * @brief Standard normal probability density function (mean 0, sd 1).
 *
 * Evaluates the standard normal PDF at z.
 *
 * @param z Standardized coordinate
 * @return PDF value
 */
double normalpdf_std(double z) {
    return exp(-pow(z, 2) / 2) / sqrt(2 * M_PI);
}

/**
 * @brief Normal distribution PDF with mean mu and standard deviation sigma.
 *
 * Evaluates the PDF at x for N(mu, sigma^2).
 *
 * @param x Point to evaluate
 * @param mu Mean
 * @param sigma Standard deviation (must be > 0)
 * @return PDF value
 */
double normalpdf(double x, double mu, double sigma) {
    double z = (x - mu) / sigma;
    return (1 / sigma) * normalpdf_std(z);
}

/**
 * @brief Standard normal cumulative distribution function (Phi).
 *
 * Returns Phi(upper) = P(Z <= upper) for Z ~ N(0,1) using erf.
 *
 * @param upper Upper limit
 * @return CDF value
 */
double Phi(double upper) {
    return 0.5 * (1 + erf(upper / sqrt(2)));
}

/**
 * @brief Standard normal probability between two z-values.
 *
 * Computes P(lower <= Z <= upper) for Z ~ N(0,1).
 *
 * @param lower Lower z-value
 * @param upper Upper z-value
 * @return Probability mass between lower and upper
 */
double normalcdf_std(double lower, double upper) {
    double below_lower = Phi(lower);
    double above_upper = 1 - Phi(upper);
    return 1 - (below_lower + above_upper);
}

/**
 * @brief Normal CDF between two values for N(mu, sigma^2).
 *
 * Standardizes inputs and returns P(lower <= X <= upper).
 *
 * @param lower Lower bound
 * @param upper Upper bound
 * @param mu Mean
 * @param sigma Standard deviation
 * @return Probability P(lower <= X <= upper)
 */
double normalcdf(double lower, double upper, double mu, double sigma) {
    double z_lower = (lower - mu) / sigma;
    double z_upper = (upper - mu) / sigma;
    return normalcdf_std(z_lower, z_upper);
}

