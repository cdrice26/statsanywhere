#include "quadrature_weights.h"

/**
 * @brief Estimate an integral using Gaussian quadrature from Jacobi matrix.
 *
 * Constructs the n-by-n Jacobi matrix from generator arrays a and b, finds
 * its eigenvalues (nodes) and eigenvectors to compute weights, and then
 * evaluates the provided integrand f at the nodes. Returns the weighted sum.
 *
 * @param n Number of quadrature nodes
 * @param a Array of length n containing diagonal generators
 * @param b Array of length n-1 containing squared off-diagonal generators
 * @param mu_0 Zeroth moment / scaling factor for weight computation
 * @param f Integrand function of a node and user context
 * @param ctx User-provided context forwarded to f
 * @return Estimated integral value, or negative on failure
 */
double quadrature_estimate(int n, double* a, double* b, double mu_0, double (*f)(double, void*), void* ctx) {
    double** J = alloc_matrix(n, n);
    if (J == NULL) return -1;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                J[i][j] = a[i];
            } else if (abs(i - j) == 1) {
                if (i < j) {
                    J[i][j] = sqrt(b[i]);
                } else {
                    J[i][j] = sqrt(b[j]);
                }
            } else {
                J[i][j] = 0;
            }
        }
    }
    EigenResult eigen_result = eigendecompose_tridiagonal(J, n, 1e-12);

    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        double first_component = eigen_result.eigenvectors[i][0];
        double weight = mu_0 * pow(first_component, 2);
        double x_i = eigen_result.eigenvalues[i];
        sum += weight * f(x_i, ctx);
    }

    free_matrix(J, n);
    free_eigen_result(&eigen_result);
    return sum;
}
