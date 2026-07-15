#include "vectors.h"
#include "matrices.h"
#include <stdbool.h>

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

double** copy_matrix(double* const* A, int m, int n) {
    double** result = alloc_matrix(m, n);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = A[i][j];
        }
    }
    return result;
}

double** transpose(double* const* A, int m, int n) {
    double** AT = alloc_matrix(m, n);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            AT[j][i] = A[i][j];
        }
    }
    return AT;
}

double** gram_schmidt(double* const* A, int n) {
    double** Q = malloc(n * sizeof(double*));
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

QR_Decomposition QR_decompose(double* const* A, int n) {
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

double* diag(double* const* A, int n) {
    double* result = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        result[i] = A[i][i];
    }
    return result;
}

double** add_matrices(double* const* A, int ma, int na, double* const* B, int mb, int nb) {
    if (ma != mb || na != nb) {
        return NULL;
    }

    double** result = alloc_matrix(ma, na);
    if (result == NULL) {
        return NULL;
    }

    for (int i = 0; i < ma; i++) {
        for (int j = 0; j < na; j++) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }

    return result;
}

double** subtract_matrices(double* const* A, int ma, int na, double* const* B, int mb, int nb) {
    if (ma != mb || na != nb) {
        return NULL;
    }

    double** result = alloc_matrix(ma, na);
    if (result == NULL) {
        return NULL;
    }

    for (int i = 0; i < ma; i++) {
        for (int j = 0; j < na; j++) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }

    return result;
}


double** multiply(double* const* A, int ma, int na, double* const* B, int mb, int nb) {
    if (na != mb) {
        return NULL;
    }

    double** result = alloc_matrix(ma, nb);
    if (result == NULL) {
        return NULL;
    }

    for (int i = 0; i < ma; i++) {
        for (int j = 0; j < nb; j++) {
            double sum = 0.0;
            for (int k = 0; k < na; k++) {
                sum += A[i][k] * B[k][j];
            }
            result[i][j] = sum;
        }
    }

    return result;
}

double* eigenvalues(double* const* A, int n, double tol) {
    double** A_copy = copy_matrix(A, n, n);

    for (int i = 0; i < 20; i++) {
        QR_Decomposition qr = QR_decompose(A_copy, n);

        double** newA = multiply(qr.R, n, n, qr.Q, n, n);
        free_matrix(qr.Q, n);
        free_matrix(qr.R, n);

        double** diff = subtract_matrices(newA, n, n, A, n, n);
        double* diagonal = diag(diff, n);
        free_matrix(diff, n);

        for (int j = 0; j < n; j++) {
            diagonal[j] = fabs(diagonal[j]);
        }

        double m = max(diagonal, n);
        free(diagonal);

        if (m < tol) {
            free_matrix(A_copy, n);
            A_copy = newA;
            break;
        }

        free_matrix(A_copy, n);
        A_copy = newA;
    }

    double* eigenvalues = diag(A_copy, n);
    free_matrix(A_copy, n);
    return eigenvalues;
}
