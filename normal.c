#include "normal.h"

typedef struct {
    double mu;
    double sigma;
    double (*user_f)(double);
} NormalContext;

static double normal_transform_and_call(double t, void* ctx) {
    NormalContext* nc = (NormalContext*)ctx;
    double x = nc->mu + nc->sigma * t;
    return nc->user_f(x);
}

double normal_expectation(double mu, double sigma, int n, double (*f)(double)) {
    double* a = malloc(n * sizeof(double));
    double* b = malloc((n - 1) * sizeof(double));

    for (int i = 0; i < n; i++) a[i] = 0.0;
    for (int i = 0; i < n - 1; i++) b[i] = (double)(i + 1);

    double beta_0 = sqrt(2.0 * M_PI);
    NormalContext ctx = { mu, sigma, f };

    double raw = quadrature_estimate(n, a, b, beta_0, normal_transform_and_call, &ctx);

    free(a);
    free(b);

    return raw / beta_0;
}

double normalpdf_std(double z) {
    return exp(-pow(z, 2) / 2) / sqrt(2 * M_PI);
}

double normalpdf(double x, double mu, double sigma) {
    double z = (x - mu) / sigma;
    return (1 / sigma) * normalpdf_std(z);
}

double normalcdf(double lower, double upper, double mu, double sigma) {
    double multiplier = 1.0 / sqrt(M_PI);
    return 0.0;
    /* calculate J at several n's, check for convergence. n = 5
    prev = gauss_quadrature_estimate(n)
    loop:
        n = n * 2
        curr = gauss_quadrature_estimate(n)
        if abs(curr - prev) < tol:
            return curr
        prev = curr

        hermite, so a = 0 and b = n. */
}

double one(double x) { return 1.0; }
double ident(double x) { return x; }
double square(double x) { return x * x; }

int main() {
    printf("%f", normal_expectation(5.0, 2.0, 20, one));     //should return ~1.0        (total probability)
    printf("%f", normal_expectation(5.0, 2.0, 20, ident));   //should return ~5.0        (the mean, mu)
    printf("%f", normal_expectation(5.0, 2.0, 20, square));  //should return ~29.0       (E[X^2] = sigma^2 + mu^2 = 4+25)
}
