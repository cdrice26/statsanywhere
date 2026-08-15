#include "binomial.h"

/**
 * @brief Binomial probability mass function.
 *
 * Computes P(X = k) for X ~ Binomial(n, p).
 *
 * @param n Number of trials
 * @param k Number of successes
 * @param p Probability of success
 * @return Probability mass P(X = k)
 */
double binompmf(int n, int k, double p) {
    return combination(n, k) * pow(p, k) * pow(1 - p, n - k);
}

/**
 * @brief Binomial cumulative distribution function (lower tail).
 *
 * Sums binompdf from 0 to floor(k) to compute P(X <= k) for
 * X ~ Binomial(n, p).
 *
 * @param n Number of trials
 * @param k Threshold (may be non-integer; floor(k) is used)
 * @param p Probability of success
 * @return Cumulative probability P(X <= k)
 */
double binomcdf(int n, double k, double p) {
    double result = 0;
    for (int i = 0; i <= floor(k); i++) {
        result += binompmf(n, i, p);
    }
    return result;
}

int main() {
    printf("binompdf(5, 3, .7) (should be 0.308700) = %f\n", binompmf(5, 3, 0.7));
    printf("binomcdf(6, 2, .7) (should be 0.070470) = %f\n", binomcdf(6, 2, 0.7));
    printf("binompdf(100, 50, .7) (should be 0.000013) = %f\n", binompmf(100, 50, 0.7));
    printf("binomcdf(100, 50, .7) (should be 0.000022) = %f\n", binomcdf(100, 50, 0.7));
}
