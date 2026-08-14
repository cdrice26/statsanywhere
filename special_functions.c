#include "special_functions.h"
#include "quadrature_weights.h"

typedef struct {
    double z;
} GammaCtx;

typedef struct {
    double s;
    double x;
} LowerIncompleteGammaCtx;

typedef struct {
    double x;
    double a;
    double b;
} IncompleteBetaCtx;

double gamma_integrand(double t, void* ctx) {
    double z = ((GammaCtx*)ctx)->z;
    return pow(t, z - 1.0);
}

double lower_incomplete_gamma_integrand(double t, void* ctx) {
    LowerIncompleteGammaCtx* gctx = (LowerIncompleteGammaCtx*)ctx;
    double s = gctx->s;
    double x = gctx->x;
    return exp(-(x / 2) * (t + 1));
}

double incomplete_beta_integrand(double t, void* ctx) {
    IncompleteBetaCtx* bctx = (IncompleteBetaCtx*)ctx;
    double x = bctx->x;
    double b = bctx->b;
    return pow(1 - (x / 2.0) * (t + 1.0), b - 1);
}

void set_jacobi_matrix_generators_incomplete(double* alpha, double* beta_sq, int n, double s) {
    for (int i = 0; i < n; i++) {
        double v = i + 1;
        if (fabs(s - 1) < 1e-12) {
            alpha[i] = 0.0;
        } else {
            alpha[i] = pow(s - 1, 2) / ((2 * v + s - 3) * (2 * v + s - 1));
        }
        if (i != n - 1) {
            beta_sq[i] = (4 * pow(v, 2) * pow(v + s - 1, 2)) /
                (pow(2 * v + s - 1, 2) * (2 * v + s - 2) * (2 * v + s));
        }
    }
}

double gamma(double z) {
    GammaCtx ctx = { z };
    int n = 50;
    double* alpha = malloc(n * sizeof(double));
    double* beta_sq = malloc((n - 1) * sizeof(double));
    for (int i = 0; i < n; i++) {
        alpha[i] = 2.0 * (i + 1) - 1.0;
        if (i != n - 1) {
            beta_sq[i] = pow(i + 1, 2);
        }
    }
    double integral = quadrature_estimate(n, alpha, beta_sq, 1.0, gamma_integrand, &ctx);
    free(alpha);
    free(beta_sq);
    return integral;
}

double lower_incomplete_gamma(double s, double x) {
    LowerIncompleteGammaCtx ctx = { s, x };
    int n = 50;
    double* alpha = malloc(n * sizeof(double));
    double* beta_sq = malloc((n - 1) * sizeof(double));
    set_jacobi_matrix_generators_incomplete(alpha, beta_sq, n, s);
    double outside_factor = pow(x / 2, s);
    double mu_0 = pow(2, s) / s;
    double integral = quadrature_estimate(n, alpha, beta_sq, mu_0, lower_incomplete_gamma_integrand, &ctx);
    free(alpha);
    free(beta_sq);
    return integral * outside_factor;
}

double incomplete_beta(double x, double a, double b) {
    IncompleteBetaCtx ctx = { x, a, b };
    int n = 50;
    double* alpha = malloc(n * sizeof(double));
    double* beta_sq = malloc((n - 1) * sizeof(double));
    set_jacobi_matrix_generators_incomplete(alpha, beta_sq, n, a);
    double outside_factor = pow(x / 2.0, a);
    double mu_0 = pow(2, a) / a;
    double integral = quadrature_estimate(n, alpha, beta_sq, mu_0, incomplete_beta_integrand, &ctx);
    free(alpha);
    free(beta_sq);
    return integral * outside_factor;
}

double regularized_incomplete_beta(double x, double a, double b) {
    return incomplete_beta(x, a, b) / beta(a, b);
}

double beta(double z1, double z2) {
    return (gamma(z1) * gamma(z2)) / gamma(z1 + z2);
}
