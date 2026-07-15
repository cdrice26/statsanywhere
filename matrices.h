#ifndef MATRICES_H
#define MATRICES_H
#include <stdlib.h>

typedef struct {
    double** Q;
    double** R;
} QR_Decomposition;

double** alloc_matrix(int m, int n);
void free_matrix(double** A, int m);
double** copy_matrix(double* const* A, int m, int n);
double** transpose(double* const* A, int m, int n);
double** gram_schmidt(double* const* A, int n);
QR_Decomposition QR_decompose(double* const* A, int n);
double* diag(double* const* A, int n);
double** add_matrices(double* const* A, int ma, int na, double* const* B, int mb, int nb);
double** subtract_matrices(double* const* A, int ma, int na, double* const* B, int mb, int nb);
double** multiply(double* const* A, int ma, int na, double* const* B, int mb, int nb);
double max(const double* v, int n);
double* eigenvalues(double* const* A, int n, double tol);

#endif
