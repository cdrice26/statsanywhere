#include "quadrature_weights.h"

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
