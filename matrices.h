#ifndef MATRICES_H
#define MATRICES_H
#include <stdlib.h>

typedef struct {
    double** Q;
    double** R;
} QR_Decomposition;

double** alloc_matrix(int m, int n);
void free_matrix(double** A, int m);
double** copy_matrix(double** A, int m, int n);
double** transpose(double** A, int m, int n);
double** gram_schmidt(double** A, int n);
QR_Decomposition QR_decompose(double** A, int n);

#endif
