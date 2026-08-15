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

