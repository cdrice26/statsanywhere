#include "vectors.h"
#include "matrices.h"

double** alloc_matrix(int m, int n) {
    double** A = (double**)malloc(m * sizeof(double*));
    for (int i = 0; i < m; i++) {
        A[i] = (double*)malloc(n * sizeof(double));
    }
    return A;
}

void free_matrix(double** A, int m) {
    for (int i = 0; i < m; i++) {
        free(A[i]);
    }
    free(A);
}

double** copy_matrix(double** A, int m, int n) {
    double** result = alloc_matrix(m, n);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = A[i][j];
        }
    }
    return result;
}

double** transpose(double** A, int m, int n) {
    double** AT = alloc_matrix(m, n);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            AT[j][i] = A[i][j];
        }
    }
    return AT;
}

double** gram_schmidt(double** A, int n) {
    double** Q = alloc_matrix(n, n);
    for (int i = 0; i < n; i++) {
        Q[i] = copy_vector(A[i], n);
        for (int j = 0; j < i; j++) {
            double* proj = projection(A[i], Q[j], n);
            double* newQ = subtract(Q[i], proj, n);

            free(proj);
            free(Q[i]);
            Q[i] = newQ;
        }
    }
    for (int i = 0; i < n; i++) {
        double mag = magnitude(Q[i], n);
        double* newQ = multiply_scalar(Q[i], n, 1.0 / mag);

        free(Q[i]);
        Q[i] = newQ;
    }
    return Q;
}

QR_Decomposition QR_decompose(double** A, int n) {
    double** AT = transpose(A, n, n);
    double** Qtmp = gram_schmidt(AT, n);
    double** Q = transpose(Qtmp, n, n);
    double** R = alloc_matrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i > j) {
                R[i][j] = 0.0;
            } else {
                R[i][j] = dot_product(Qtmp[i], AT[j], n);
            }
        }
    }
    free_matrix(AT, n);
    free_matrix(Qtmp, n);
    QR_Decomposition result = {Q, R};
    return result;
}
