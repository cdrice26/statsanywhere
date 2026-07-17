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

double** multiply_matrix_by_scalar(double* const* A, int rows, int cols, double s) {
    double** R = malloc(rows * sizeof *R);
    if (!R) return NULL;

    for (int i = 0; i < rows; i++) {
        R[i] = malloc(cols * sizeof *R[i]);
        if (!R[i]) {
            for (int k = 0; k < i; k++) free(R[k]);
            free(R);
            return NULL;
        }
        for (int j = 0; j < cols; j++) {
            R[i][j] = A[i][j] * s;
        }
    }
    return R;
}

double** identity(int n) {
    double** A = malloc(n * sizeof *A);
    if (!A) return NULL;

    for (int i = 0; i < n; i++) {
        A[i] = malloc(n * sizeof *A[i]);
        if (!A[i]) {
            for (int k = 0; k < i; k++) free(A[k]);
            free(A);
            return NULL;
        }
        for (int j = 0; j < n; j++) {
            A[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
    return A;
}

void swap_rows(double **A, int i, int j) {
    if (i == j) return;
    double *tmp = A[i];
    A[i] = A[j];
    A[j] = tmp;
}

void scale_row(double **A, int cols, int i, double c) {
    for (int k = 0; k < cols; ++k) A[i][k] *= c;
}

void add_multiple_row(double **A, int cols, int i, int j, double c) {
    if (c == 0.0) return;
    for (int k = 0; k < cols; ++k) A[i][k] += c * A[j][k];
}

double** reduce_back(EliminationResult* result) {
    double** A = copy_matrix(result->rref, result->m, result->n);
    if (A == NULL) return NULL;
    Node* pivot_ptr = result->pivots;
    while (pivot_ptr->next != NULL) {
        for (int i = pivot_ptr->row - 1; i >= 0; i--) {
            add_multiple_row(A, result->n, i, pivot_ptr->row, -A[i][pivot_ptr->col]);
        }
        pivot_ptr = pivot_ptr->next;
    }
    return A;
}

EliminationResult reduce(double* const* A, int m, int n, double tol) {
    Node* pivot_list = create_linked_list();
    double** B = copy_matrix(A, m, n);

    if (pivot_list == NULL || B == NULL) {
        free_linked_list(pivot_list);
        EliminationResult result = { NULL, 0, 0, NULL };
        return result;
    }

    // Step 5: Repeat the forward phase for all rows of the matrix by looping
    // through the rows of A
    for (int i = 0; i < m; i++) {
        int col1 = -1;
        int row1 = -1;

        // Step 1: Find the leftmost column that does not consist entirely of zeros.
        // Call this column 'col1'. If none exists, end the algorithm.
        for (int j = i; j < n; j++) {
            bool should_break = false;
            for (int k = i; k < m; k++) {
                if (fabs(B[k][j]) > tol) {
                    col1 = j;
                    should_break = true;
                    break;
                }
            }
            if (should_break) break;
        }
        if (col1 == -1) break;

        // Step 2: Bring a non-zero entry to the top of col1
        for (int j = i; j < m; j++) {
            if (fabs(B[j][col1]) > tol) {
                row1 = j;
                break;
            }
        }
        if (row1 == -1) {
            break;
        }
        swap_rows(B, i, row1);

        // Step 3: Make the entry in the top of col1 1
        if (fabs(B[i][col1]) > tol) {
            scale_row(B, n, i, 1.0 / B[i][col1]);
        }
        push(pivot_list, i, col1);

        // Step 4: Make all entries below A[i,col1] 0's
        for (int j = i + 1; j < m; j++) {
            if (fabs(B[j][col1]) > tol) {
                add_multiple_row(B, n, j, i, -B[j][col1]);
            }
        }
    }

    EliminationResult ref = {B, m, n, pivot_list};
    EliminationResult result = {reduce_back(&ref), m, n, pivot_list};
    free_matrix(B, m);
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

double** eigenspace(double* const* A, int n, double tol) {
    double* lambdas = eigenvalues(A, n, tol);
}
