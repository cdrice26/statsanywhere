#include "special_functions.h"
#include "quadrature_weights.h"

typedef struct {
    double z;
} GammaCtx;

typedef struct {
    double s;
    double x;
} LowerIncompleteGammaCtx;

// Original code using Gauss-Legendre - not precise enough

// double gamma_integrand(double x, void* ctx) {
//     double z = ((GammaCtx*)ctx)->z;
//     return pow(-log((x + 1.0) / 2.0), z - 1.0);
// }

// double gamma(double z) {
//     GammaCtx ctx = { z };
//     int n = 300;
//     double* alpha = malloc(n * sizeof(double));
//     double* beta_sq = malloc((n - 1) * sizeof(double));
//     for (int i = 0; i < n; i++) {
//         alpha[i] = 0;
//         if (i != n - 1) {
//             beta_sq[i] = pow(i + 1, 2) / (4 * pow(i + 1, 2) - 1);
//         }
//     }
//     double integral = quadrature_estimate(n, alpha, beta_sq, 2.0, gamma_integrand, &ctx);
//     return 0.5 * integral;
// }

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
    double outside_factor = pow(x / 2, s);
    double mu_0 = pow(2, s) / s;
    double integral = quadrature_estimate(n, alpha, beta_sq, mu_0, lower_incomplete_gamma_integrand, &ctx);
    free(alpha);
    free(beta_sq);
    return integral * outside_factor;
}

int main() {
    printf("(Should be 24) Gamma(5) = %f\n", gamma(5.0));
    printf("(Should be 413.40751) Gamma(6.7) = %f\n", gamma(6.7));
    printf("(Should be 0.9933) LGamma(1, 5) = %f\n", lower_incomplete_gamma(1, 5));
    printf("(Should be 1.8232) LGamma(3, 5.5) = %f\n", lower_incomplete_gamma(3, 5.5));
}
