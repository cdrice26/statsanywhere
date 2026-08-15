#include "operations.h"

/**
 * @brief Compute the factorial of n.
 *
 * Returns n! as a double. Uses iterative multiplication. Note: the
 * implementation does not perform overflow checks for large n.
 *
 * @param n Non-negative integer
 * @return n! as a double
 */
double factorial(unsigned int n) {
    unsigned int result = n;
    for (int i = n - 1; i > 0; i--) {
        result *= i;
    }
    return result;
}

/**
 * @brief Compute the binomial coefficient C(n, r).
 *
 * Calculates "n choose r" using a multiplicative approach that reduces
 * intermediate values by using the smaller of r and n-r. Returns 0 when r
 * > n. Results are returned as a double.
 *
 * @param n Total number of items
 * @param r Number of selections
 * @return Binomial coefficient C(n, r) as a double
 */
double combination(unsigned int n, unsigned int r) {
    if (r > n || r < 0) return 0;
    if (r == 0 || r == n) return 1;

    // Use the smaller of r and n-r to minimize computation
    if (r > n - r) r = n - r;

    double result = 1;
    for (unsigned int i = 0; i < r; i++) {
        result *= (double)(n - i) / (double)(i + 1);
    }
    return result;
}
