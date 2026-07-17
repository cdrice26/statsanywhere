#ifndef MATRICES_H
#define MATRICES_H
#include <stdlib.h>
#include "linked_list.h"
#include "vectors.h"
#include <stdbool.h>
#include <stdio.h>

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

typedef struct {
    double* eigenvalues;   // length n
    double** eigenvectors; // array of `count` vectors, each length n
    int n;                 // length of each eigenvalue's... eigenvector (dimension)
    int count;              // number of eigenvectors found (rows in eigenvectors)
} EigenResult;

typedef struct { int row, col; } Point;

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
double** partition(int m, double* const* A, int an, double* const* B, int bn);
EliminationResult reduce(double* const* A, int m, int n, int pivot_col_limit, double tol);
double** identity(int n);
double* eigenvalues(double* const* A, int n, double tol);
void free_eigen_result(EigenResult* r);
EigenResult eigenvectors(double* const* A, int n, double tol);

#endif
