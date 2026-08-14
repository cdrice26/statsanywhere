// main.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrices.h"
#include "operations.h"
#include "special_functions.h"

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

static void print_matrix_rect(const char* name, double** M, int m, int n) {
    printf("%s =\n", name);
    for (int i = 0; i < m; i++) {
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
    A[0][0] = 2;  A[0][1] = -1; A[0][2] = 4;
    A[1][0] = -1;   A[1][1] = 3; A[1][2] = 0;
    A[2][0] = 4;  A[2][1] = 0;  A[2][2] = 5;

    EliminationResult res = reduce(A, n, n, n, 1e-6);

    print_matrix("RREF", res.rref, n);

    QR_Decomposition qr = QR_decompose(A, n);

    print_matrix("Q", qr.Q, n);
    print_matrix("R", qr.R, n);

    double** QR = mul_square(qr.Q, qr.R, n);
    print_matrix("Q*R", QR, n);

    printf("max |(Q*R) - A| = %.12f\n", max_abs_diff(QR, A, n));

    EigenResult eigen = eigendecompose(A, n, 1e-6);
    print_array(eigen.eigenvalues, n);
    print_matrix_rect("Eigenvectors", eigen.eigenvectors, eigen.count, eigen.n);

    // Cleanup
    free_matrix(A, n);
    free_matrix(qr.Q, n);
    free_matrix(qr.R, n);
    free_matrix(QR, n);
    free_matrix(eigen.eigenvectors, eigen.count);
    free(eigen.eigenvalues);

    printf("(should be 362,880) 9! = %f\n", factorial(9));
    printf("(should be 15) 6C4 = %f\n", combination(6, 4));

    printf("(Should be 24) Gamma(5) = %f\n", gamma(5.0));
    printf("(Should be 413.40751) Gamma(6.7) = %f\n", gamma(6.7));
    printf("(Should be 0.9933) LGamma(1, 5) = %f\n", lower_incomplete_gamma(1, 5));
    printf("(Should be 1.8232) LGamma(3, 5.5) = %f\n", lower_incomplete_gamma(3, 5.5));
    printf("(Should be 0.0833) Beta(3, 2) = %f\n", beta(3, 2));
    printf("(Should be 0.00035) Beta(4, 10) = %f\n", beta(4, 10));

    return 0;
}
