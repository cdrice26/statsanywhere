#include "special_functions.h"
#include "quadrature_weights.h"

typedef struct {
    double z;
} GammaCtx;

typedef GammaCtx ErrorCtx;

typedef struct {
    double s;
    double x;
} LowerIncompleteGammaCtx;

typedef struct {
    double x;
    double a;
    double b;
} IncompleteBetaCtx;

/**
 * @brief Integrand for the Gamma function quadrature.
 *
 * Returns t^(z-1) for use with the quadrature scheme. Context provides z.
 *
 * @param t Integration variable
 * @param ctx Pointer to GammaCtx containing z
 * @return Evaluated integrand value
 */
double gamma_integrand(double t, void* ctx) {
    double z = ((GammaCtx*)ctx)->z;
    return pow(t, z - 1.0);
}

/**
 * @brief Integrand for the lower incomplete gamma quadrature.
 *
 * Uses parameters s and x from context to evaluate exp(-(x/2)*(t+1)).
 *
 * @param t Integration variable
 * @param ctx Pointer to LowerIncompleteGammaCtx containing s and x
 * @return Evaluated integrand value
 */
double lower_incomplete_gamma_integrand(double t, void* ctx) {
    LowerIncompleteGammaCtx* gctx = (LowerIncompleteGammaCtx*)ctx;
    double s = gctx->s;
    double x = gctx->x;
    return exp(-(x / 2) * (t + 1));
}

/**
 * @brief Integrand for the incomplete beta quadrature.
 *
 * Uses a context with x and b to compute (1 - (x/2)*(t+1))^(b-1).
 *
 * @param t Integration variable
 * @param ctx Pointer to IncompleteBetaCtx containing x, a, b
 * @return Evaluated integrand value
 */
double incomplete_beta_integrand(double t, void* ctx) {
    IncompleteBetaCtx* bctx = (IncompleteBetaCtx*)ctx;
    double x = bctx->x;
    double b = bctx->b;
    return pow(1 - (x / 2.0) * (t + 1.0), b - 1);
}

/**
 * @brief Integrand for the error function quadrature.
 *
 * Computes exp(-((z/2)*(t+1))^2) using z from context.
 *
 * @param t Integration variable
 * @param ctx Pointer to ErrorCtx containing z
 * @return Evaluated integrand value
 */
double erf_integrand(double t, void* ctx) {
    ErrorCtx* ectx = (ErrorCtx*)ctx;
    double z = ectx->z;
    return exp(-pow((z / 2.0) * (t + 1), 2));
}

/**
 * @brief Populate Jacobi matrix generator arrays for incomplete integrals.
 *
 * Fills alpha and beta_sq for n nodes using parameter s.
 *
 * @param alpha Output array of length n for diagonal generators
 * @param beta_sq Output array of length n-1 for squared off-diagonals
 * @param n Number of nodes
 * @param s Parameter (e.g., shape) controlling generator values
 */
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

/**
 * @brief Compute the Gamma function via quadrature.
 *
 * Uses Gauss-Laguerre quadrature with Jacobi matrix eigenvalues/weights to
 * estimate the Gamma function integral. Allocates temporary arrays and
 * returns the integral value.
 *
 * @param z Argument to the Gamma function
 * @return Approximated Gamma(z)
 */
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

/**
 * @brief Compute the lower incomplete gamma function.
 *
 * Uses the Gauss-Jacobi quadrature with Jacobi generators tailored to parameter s
 * to approximate the lower incomplete gamma. Allocates temporaries and
 * returns the evaluated value.
 *
 * @param s Shape parameter
 * @param x Integration upper limit parameter
 * @return Approximated lower incomplete gamma value
 */
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

/**
 * @brief Compute the incomplete beta function.
 *
 * Uses Gauss-Jacobi quadrature with Jacobi generators derived from parameter a to
 * estimate the incomplete beta integral. Returns the evaluated value.
 *
 * @param x Upper integration limit parameter
 * @param a First shape parameter
 * @param b Second shape parameter
 * @return Approximated incomplete beta value (unnormalized)
 */
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

/**
 * @brief Compute the regularized incomplete beta function.
 *
 * Returns the incomplete beta divided by the complete Beta(a,b).
 *
 * @param x Upper integration limit parameter
 * @param a First shape parameter
 * @param b Second shape parameter
 * @return Regularized incomplete beta in [0,1]
 */
double regularized_incomplete_beta(double x, double a, double b) {
    return incomplete_beta(x, a, b) / beta(a, b);
}

/**
 * @brief Compute the Beta function via Gamma functions.
 *
 * Uses the identity Beta(z1,z2) = Gamma(z1)*Gamma(z2)/Gamma(z1+z2).
 *
 * @param z1 First argument
 * @param z2 Second argument
 * @return Beta(z1, z2)
 */
double beta(double z1, double z2) {
    return (gamma(z1) * gamma(z2)) / gamma(z1 + z2);
}

/**
 * @brief Compute the error function erf(z).
 *
 * Estimates erf using Gauss-Legendre quadrature with precomputed Jacobi generators and
 * returns the scaled integral value.
 *
 * @param z Argument to erf
 * @return Approximated erf(z)
 */
double erf(double z) {
    ErrorCtx ctx = { z };
    int n = 50;
    double* alpha = malloc(n * sizeof(double));
    double* beta_sq = malloc((n - 1) * sizeof(double));
    for (int i = 0; i < n; i++) {
        alpha[i] = 0;
        if (i != n - 1) {
            beta_sq[i] = pow(i + 1, 2) / (4 * pow(i + 1, 2) - 1);
        }
    }
    double outside_factor = z / sqrt(M_PI);
    double integral = quadrature_estimate(n, alpha, beta_sq, 2.0, erf_integrand, &ctx);
    free(alpha);
    free(beta_sq);
    return integral * outside_factor;
}
