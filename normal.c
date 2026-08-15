#include "normal.h"

typedef struct {
    double max;
} MaxVal;

double normalpdf_std(double z) {
    return exp(-pow(z, 2) / 2) / sqrt(2 * M_PI);
}

double normalpdf(double x, double mu, double sigma) {
    double z = (x - mu) / sigma;
    return (1 / sigma) * normalpdf_std(z);
}

double Phi(double upper) {
    return 0.5 * (1 + erf(upper / sqrt(2)));
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
