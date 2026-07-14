#include "vectors.h"

double* add(double* v1, double* v2, int n) {
    double* result = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        result[i] = v1[i] + v2[i];
    }
    return result;
}

double* negate(double* v, int n) {
    double* result = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        result[i] = -v[i];
    }
    return result;
}

double* subtract(double* v1, double* v2, int n) {
    double* nv2 = negate(v2, n);
    double* result = add(v1, nv2, n);
    free(nv2);
    return result;
}

double dot_product(double* v1, double* v2, int n) {
    double result = 0;
    for (int i = 0; i < n; i++) {
        result += v1[i] * v2[i];
    }
    return result;
}

double* multiply_scalar(double* v, int n, double k) {
    double* result = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        result[i] = v[i] * k;
    }
    return result;
}

double* projection(double* of, double* onto, int n) {
    return multiply_scalar(onto, n, dot_product(onto, of, n) / dot_product(onto, onto, n));
}

double magnitude(double* v, int n) {
    return sqrt(dot_product(v, v, n));
}

double* copy_vector(double* v, int n) {
    return multiply_scalar(v, n, 1.0);
}
