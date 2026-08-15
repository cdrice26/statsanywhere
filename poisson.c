#include "poisson.h"

/**
 * @brief Poisson probability mass function.
 *
 * Computes P(K = k) for K ~ Poisson(lambda).
 *
 * @param k Number of occurrences (integer)
 * @param lambda Rate parameter (>= 0)
 * @return Probability mass P(K = k)
 */
double poissonpmf(int k, double lambda) {
    return (pow(lambda, k) * exp(-lambda)) / factorial(k);
}

/**
 * @brief Poisson cumulative distribution function (lower tail).
 *
 * Computes P(K <= k) by summing terms of the Poisson PMF up to floor(k).
 * Uses an incremental log-factorial to stabilize large-k computations.
 *
 * @param k Threshold (may be non-integer; floor(k) is used)
 * @param lambda Rate parameter
 * @return Cumulative probability P(K <= k)
 */
double poissoncdf(double k, double lambda) {
    double result = 0;
    double log_factorial_term = 0.0;

    for (int j = 0; j <= floor(k); j++) {
        double term = exp(j * log(lambda) - log_factorial_term);
        result += term;
        log_factorial_term += log(j + 1);
    }
    return result * exp(-lambda);
}

int main() {
    printf("poissonpdf(3, 5) (should be 0.140374) = %f\n", poissonpmf(3, 5));
    printf("poissoncdf(3, 5) (should be 0.265026) = %f\n", poissoncdf(3, 5));
}
