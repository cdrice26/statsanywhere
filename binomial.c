#include "binomial.h"

double binompdf(int n, int k, double p) {
    return combination(n, k) * pow(p, k) * pow(1 - p, n - k);
}

double binomcdf(int n, double k, double p) {
    double result = 0;
    for (int i = 0; i <= floor(k); i++) {
        result += binompdf(n, i, p);
    }
    return result;
}

int main() {
    printf("binompdf(5, 3, .7) (should be 0.308700) = %f\n", binompdf(5, 3, 0.7));
    printf("binomcdf(6, 2, .7) (should be 0.070470) = %f", binomcdf(6, 2, 0.7));
}
