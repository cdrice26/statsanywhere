#include "poisson.h"

double poissonpdf(int k, double lambda) {
    return (pow(lambda, k) * exp(-lambda)) / factorial(k);
}

double poissoncdf(double k, double lambda) {
    double result = 0;
    double log_factorial_term = 0.0;

    for (int j = 0; j <= floor(k); j++) {
        double term = exp(j * log(lambda) - log_factorial_term);
        result += term;
        log_factorial_term += log(j + 1);
    }
    return result * exp(-lambda);
}

int main() {
    printf("poissonpdf(3, 5) (should be 0.140374) = %f\n", poissonpdf(3, 5));
    printf("poissoncdf(3, 5) (should be 0.265026) = %f\n", poissoncdf(3, 5));
}
