#ifndef MATRICES_H
#define MATRICES_H
#include <stdlib.h>
#include "linked_list.h"
#include "vectors.h"
#include <stdbool.h>

typedef struct {
    double** Q;
    double** R;
} QR_Decomposition;

typedef struct {
    double** rref;
    int m;
    int n;
    Node* pivots;
} EliminationResult;

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
double** multiply_matrix_by_scalar(double* const* A, int rows, int cols, double s);
double max(const double* v, int n);
EliminationResult reduce(double* const* A, int m, int n, double tol);
double** identity(int n);
double* eigenvalues(double* const* A, int n, double tol);

#endif
