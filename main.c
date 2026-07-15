// main.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrices.h"

static void print_array(const double a[], int n) {
    for (int i = 0; i < n; i++) {
        printf("%f%s", a[i], (i + 1 == n) ? "" : " ");
    }
    printf("\n");
}

static void print_matrix(const char* name, double** M, int n) {
    printf("%s =\n", name);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("% .10f ", M[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

static double** alloc_square(int n) {
    double** M = alloc_matrix(n, n);
    return M;
}

static double** mul_square(double** A, double** B, int n) {
    double** C = alloc_matrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = sum;
        }
    }
    return C;
}

static double max_abs_diff(double** A, double** B, int n) {
    double maxv = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double d = fabs(A[i][j] - B[i][j]);
            if (d > maxv) maxv = d;
        }
    }
    return maxv;
}

int main(void) {
    int n = 3;

    // A = [ 12  -51   4
    //       6  167 -68
    //      -4   24 -41 ]
    double** A = alloc_square(n);
    A[0][0] = 12;  A[0][1] = -51; A[0][2] = 4;
    A[1][0] = 6;   A[1][1] = 167; A[1][2] = -68;
    A[2][0] = -4;  A[2][1] = 24;  A[2][2] = -41;

    QR_Decomposition qr = QR_decompose(A, n);

    print_matrix("Q", qr.Q, n);
    print_matrix("R", qr.R, n);

    double** QR = mul_square(qr.Q, qr.R, n);
    print_matrix("Q*R", QR, n);

    printf("max |(Q*R) - A| = %.12f\n", max_abs_diff(QR, A, n));

    double* eigen = eigenvalues(A, n, 1e-6);
    print_array(eigen, n);

    // Cleanup
    free(eigen);
    free_matrix(A, n);
    free_matrix(qr.Q, n);
    free_matrix(qr.R, n);
    free_matrix(QR, n);

    return 0;
}
