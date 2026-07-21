#include "vectors.h"

double* add(const double* v1, const double* v2, int n) {
    double* result = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        result[i] = v1[i] + v2[i];
    }
    return result;
}

double* negate(const double* v, int n) {
    double* result = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        result[i] = -v[i];
    }
    return result;
}

double* subtract(const double* v1, const double* v2, int n) {
    double* nv2 = negate(v2, n);
    double* result = add(v1, nv2, n);
    free(nv2);
    return result;
}

double dot_product(const double* v1, const double* v2, int n) {
    double result = 0;
    for (int i = 0; i < n; i++) {
        result += v1[i] * v2[i];
    }
    return result;
}

double* multiply_scalar(const double* v, int n, double k) {
    double* result = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        result[i] = v[i] * k;
    }
    return result;
}

double* projection(const double* of, const double* onto, int n) {
    return multiply_scalar(onto, n, dot_product(onto, of, n) / dot_product(onto, onto, n));
}

double magnitude(const double* v, int n) {
    return sqrt(dot_product(v, v, n));
}

double* copy_vector(const double* v, int n) {
    return multiply_scalar(v, n, 1.0);
}

double max(const double* v, int n) {
    double max = v[0];
    for (int i = 1; i < n; i++) {
        if (v[i] > max) max = v[i];
    }
    return max;
}

double norm(const double* v, int n) {
    double sum = 0;
    for (int i = 0; i < n; i++) {
        sum += pow(v[i], 2);
    }
    return sqrt(sum);
}
