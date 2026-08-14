#include "f.h"
#include "special_functions.h"

double fpdf(double x, double d1, double d2) {
    return (pow(d1 / d2, d1 / 2) * pow(x, d1 / 2 - 1) * pow(1 + d1 / d2 * x, -(d1 + d2) / 2)) /
    beta(d1 / 2, d2 / 2);
}

double fcdf(double x, double d1, double d2) {
    double rbeta_upper_limit = (d1 * x) / (d1 * x + d2);
    return regularized_incomplete_beta(rbeta_upper_limit, d1 / 2.0, d2 / 2.0);
}

int main() {
    printf("(Should be .071589) fpdf(3, 5, 4) = %f\n", fpdf(3, 5, 4));
    printf("(Should be .8453) fcdf(3, 5, 4) = %f\n", fcdf(3, 5, 4));
}
