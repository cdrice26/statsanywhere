#include "operations.h"

unsigned int factorial(unsigned int n) {
    unsigned int result = n;
    for (int i = n - 1; i > 0; i--) {
        result *= i;
    }
    return result;
}

double combination(unsigned int n, unsigned int r) {
    if (r > n || r < 0) return 0;
        if (r == 0 || r == n) return 1;

        // Use the smaller of r and n-r to minimize computation
        if (r > n - r) r = n - r;

        double result = 1;
        for (unsigned int i = 0; i < r; i++) {
            result *= (double)(n - i) / (double)(i + 1);
        }
        return result;
}
