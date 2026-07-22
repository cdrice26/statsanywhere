#include "normal.h"

typedef struct {
    double max;
} MaxVal;

static double transform(double t, void* ctx) {
    MaxVal* mv = (MaxVal*)ctx;
    double x = mv->max - exp(-t);
    double g = normalpdf_std(x) * exp(-t);
    return g * exp(t * t / 2.0);
}

double normal_estimate(int n, MaxVal* mv) {
    double* a = malloc(n * sizeof(double));
    double* b = malloc((n - 1) * sizeof(double));

    for (int i = 0; i < n; i++) a[i] = 0.0;
    for (int i = 0; i < n - 1; i++) b[i] = (double)(i + 1);

    double beta_0 = sqrt(2.0 * M_PI);

    double result = quadrature_estimate(n, a, b, beta_0, transform, mv);

    free(a);
    free(b);

    return result;
}

double normalpdf_std(double z) {
    return exp(-pow(z, 2) / 2) / sqrt(2 * M_PI);
}

double normalpdf(double x, double mu, double sigma) {
    double z = (x - mu) / sigma;
    return (1 / sigma) * normalpdf_std(z);
}

double Phi(double upper) {
    int n = 40;
    MaxVal mv = {upper};
    return normal_estimate(n, &mv);
}

// double normalcdf(double lower, double upper, double mu, double sigma) {
//     double multiplier = 1.0 / sqrt(M_PI);
//     return 0.0;
//     /* calculate J at several n's, check for convergence. n = 5
//     prev = gauss_quadrature_estimate(n)
//     loop:
//         n = n * 2
//         curr = gauss_quadrature_estimate(n)
//         if abs(curr - prev) < tol:
//             return curr
//         prev = curr

//         hermite, so a = 0 and b = n. */
// }

int main() {
    printf("Phi(1) (should be .841) = %f", Phi(1));
}
