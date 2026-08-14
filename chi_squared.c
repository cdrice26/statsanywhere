#include "chi_squared.h"
#include "special_functions.h"

double X2pdf(double x, int k) {
    if (x > 0) {
        return (pow(x, k / 2.0 - 1.0) * exp(-x / 2.0)) / (pow(2.0, k / 2.0) * gamma(k / 2.0));
    } else {
        return 0;
    }
}

double X2cdf(double x, int k) {
    return lower_incomplete_gamma(k / 2.0, x / 2.0) / gamma(k / 2.0);
}

int main() {
    printf("(Should be .073) X2pdf(5, 3) = %f\n", X2pdf(5, 3));
    printf("(Should be .828) X2cdf(5, 3) = %f\n", X2cdf(5, 3));
}
