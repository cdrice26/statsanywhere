#include "normal.h"

typedef struct {
    double max;
} MaxVal;

static double transform(double t, void* ctx) {
    MaxVal* mv = (MaxVal*)ctx;
    if (t < -6.0 || t > 20.0) return 0.0;
    double s = sinh(t);
    double x = mv->max - exp(-s);
    double dxdt = exp(-s) * cosh(t);
    double g = normalpdf_std(x) * dxdt;
    double val = g * exp(t * t / 2.0);
    return isfinite(val) ? val : 0.0;
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
    int n = 300;
    MaxVal mv = {upper};
    return normal_estimate(n, &mv);
}

double normalcdf_std(double lower, double upper) {
    double below_lower = Phi(lower);
    double above_upper = 1 - Phi(upper);
    return 1 - (below_lower + above_upper);
}

double normalcdf(double lower, double upper, double mu, double sigma) {
    double z_lower = (lower - mu) / sigma;
    double z_upper = (upper - mu) / sigma;
    return normalcdf_std(z_lower, z_upper);
}

int main() {
    printf("Phi(1) (should be .841) = %f\n", Phi(1));
    printf("normalcdf_std(-1, 1) (should be .683) = %f\n", normalcdf_std(-1, 1));
    printf("normalcdf(-1, 1, 1, 2) (should be .341) = %f\n", normalcdf(-1, 1, 1, 2));
}
