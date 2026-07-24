#include "matrices.h"
#include "linked_list.h"

double** alloc_matrix(int m, int n) {
    double** A = (double**)malloc(m * sizeof(double*));
    if (A == NULL) return NULL;
    for (int i = 0; i < m; i++) {
        A[i] = (double*)malloc(n * sizeof(double));
        if (A[i] == NULL) {
            for (int j = 0; j < i; j++) free(A[j]);
            free(A);
            return NULL;
        }
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
    if (result == NULL) return NULL;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = A[i][j];
        }
    }
    return result;
}

double** transpose(double* const* A, int m, int n) {
    double** AT = alloc_matrix(n, m);
    if (AT == NULL) return NULL;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            AT[j][i] = A[i][j];
        }
    }
    return AT;
}

double** gram_schmidt(double* const* A, int n) {
    double** Q = malloc(n * sizeof(double*));
    if (Q == NULL) return NULL;
    
    for (int i = 0; i < n; i++) {
        Q[i] = copy_vector(A[i], n);
        if (Q[i] == NULL) {
            for (int k = 0; k < i; k++) free(Q[k]);
            free(Q);
            return NULL;
        }
        for (int j = 0; j < i; j++) {
            double* proj = projection(A[i], Q[j], n);
            if (proj == NULL) {
                for (int k = 0; k <= i; k++) free(Q[k]);
                free(Q);
                return NULL;
            }
            double* newQ = subtract(Q[i], proj, n);
            free(proj);
            if (newQ == NULL) {
                for (int k = 0; k <= i; k++) free(Q[k]);
                free(Q);
                return NULL;
            }
            free(Q[i]);
            Q[i] = newQ;
        }
    }
    for (int i = 0; i < n; i++) {
        double mag = magnitude(Q[i], n);
        double* newQ = multiply_scalar(Q[i], n, 1.0 / mag);
        if (newQ == NULL) {
            for (int k = 0; k < n; k++) free(Q[k]);
            free(Q);
            return NULL;
        }
        free(Q[i]);
        Q[i] = newQ;
    }
    return Q;
}

QR_Decomposition QR_decompose(double* const* A, int n) {
    QR_Decomposition failed = {NULL, NULL};
    
    double** AT = transpose(A, n, n);
    if (AT == NULL) return failed;
    
    double** Qtmp = gram_schmidt(AT, n);
    if (Qtmp == NULL) {
        free_matrix(AT, n);
        return failed;
    }
    
    double** Q = transpose(Qtmp, n, n);
    if (Q == NULL) {
        free_matrix(AT, n);
        free_matrix(Qtmp, n);
        return failed;
    }
    
    double** R = alloc_matrix(n, n);
    if (R == NULL) {
        free_matrix(AT, n);
        free_matrix(Qtmp, n);
        free_matrix(Q, n);
        return failed;
    }
    
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
    if (result == NULL) return NULL;
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

double** partition(int m, double* const* A, int an, double* const* B, int bn) {
    double** result = alloc_matrix(m, an + bn);
    if (result == NULL) return NULL;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < an + bn; j++) {
            if (j < an) {
                result[i][j] = A[i][j];
            } else {
                result[i][j] = B[i][j - an];
            }
        }
    }
    return result;
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
    if (result == NULL || result->rref == NULL) return NULL;
    double** A = copy_matrix(result->rref, result->m, result->n);
    if (A == NULL) return NULL;
    Node* pivot_ptr = result->pivots;
    while (pivot_ptr->next != NULL) {
        Point* p = (Point*)pivot_ptr->data;
        for (int i = p->row - 1; i >= 0; i--) {
            add_multiple_row(A, result->n, i, p->row, -A[i][p->col]);
        }
        pivot_ptr = pivot_ptr->next;
    }
    return A;
}

EliminationResult reduce(double* const* A, int m, int n, int pivot_col_limit, double tol) {
    Node* pivot_list = create_linked_list();
    double** B = copy_matrix(A, m, n);

    if (pivot_list == NULL || B == NULL) {
        free_linked_list(pivot_list, free);
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
        for (int j = i; j < pivot_col_limit; j++) {
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
        Point *p = malloc(sizeof(Point));
        if (p == NULL) {
            free_linked_list(pivot_list, free);
            free_matrix(B, m);
            EliminationResult result = { NULL, 0, 0, NULL };
            return result;
        }
        p->row = i;
        p->col = col1;
        push(pivot_list, p);

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

double wilkinson_shift(double** A, int n) {
    double a = A[n-2][n-2];
    double b = A[n-2][n-1];
    double c = A[n-1][n-1];
    double delta = (a - c) / 2.0;
    double sign = (delta >= 0) ? 1.0 : -1.0;
    if (delta == 0.0) sign = 1.0;
    return c - sign * (b * b) / (fabs(delta) + sqrt(delta * delta + b * b));
}

EigenResult eigendecompose(double* const* A, int n, double tol) {
    EigenResult failed = {NULL, NULL, 0, 0};
    
    double** A_copy = copy_matrix(A, n, n);
    if (A_copy == NULL) return failed;
    
    double** V = identity(n);
    if (V == NULL) {
        free_matrix(A_copy, n);
        return failed;
    }

    int m = n;
    int max_iter_per_value = 1000;

    while (m > 1) {
        int converged = 0;
        for (int iter = 0; iter < max_iter_per_value; iter++) {
            double mu = wilkinson_shift(A_copy, m);

            for (int i = 0; i < m; i++) A_copy[i][i] -= mu;

            QR_Decomposition qr = QR_decompose(A_copy, m);
            if (qr.Q == NULL || qr.R == NULL) {
                free_matrix(A_copy, n);
                free_matrix(V, n);
                return failed;
            }
            
            double** RQ = multiply(qr.R, m, m, qr.Q, m, m);
            if (RQ == NULL) {
                free_matrix(qr.Q, m);
                free_matrix(qr.R, m);
                free_matrix(A_copy, n);
                free_matrix(V, n);
                return failed;
            }

            for (int i = 0; i < m; i++) RQ[i][i] += mu;

            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++)
                    A_copy[i][j] = RQ[i][j];
            free_matrix(RQ, m);

            double** Vblock = multiply(V, n, m, qr.Q, m, m);
            if (Vblock == NULL) {
                free_matrix(qr.Q, m);
                free_matrix(qr.R, m);
                free_matrix(A_copy, n);
                free_matrix(V, n);
                return failed;
            }
            
            for (int i = 0; i < n; i++)
                for (int j = 0; j < m; j++)
                    V[i][j] = Vblock[i][j];
            free_matrix(Vblock, n);

            free_matrix(qr.Q, m);
            free_matrix(qr.R, m);

            if (fabs(A_copy[m-1][m-2]) < tol) { converged = 1; break; }
        }
        m--;
        (void)converged;
    }

    double* eigenvalues = diag(A_copy, n);
    if (eigenvalues == NULL) {
        free_matrix(A_copy, n);
        free_matrix(V, n);
        return failed;
    }
    
    free_matrix(A_copy, n);

    double** eigenvectors_rows = transpose(V, n, n);
    if (eigenvectors_rows == NULL) {
        free(eigenvalues);
        free_matrix(V, n);
        return failed;
    }
    
    free_matrix(V, n);

    return (EigenResult){
        .eigenvalues = eigenvalues,
        .eigenvectors = eigenvectors_rows,
        .n = n,
        .count = n
    };
}

static double wilkinson_shift_tridiagonal(double a, double b, double c) {
    double delta = (a - c) / 2.0;
    double sign = (delta >= 0.0) ? 1.0 : -1.0;
    if (delta == 0.0) sign = 1.0;
    double denom = fabs(delta) + sqrt(delta * delta + b * b);
    return (denom != 0.0) ? c - sign * (b * b) / denom : c;
}

// Rotate rows k, k+1 of an m-wide-active row (in place), first `width` columns.
static void givens_left(double** T, int k, double c, double s, int width) {
    for (int j = 0; j < width; j++) {
        double a = T[k][j], b = T[k + 1][j];
        T[k][j]     =  c * a + s * b;
        T[k + 1][j] = -s * a + c * b;
    }
}

// Rotate columns k, k+1, over the first `rows` rows.
static void givens_right(double** T, int k, double c, double s, int rows) {
    for (int i = 0; i < rows; i++) {
        double a = T[i][k], b = T[i][k + 1];
        T[i][k]     =  c * a + s * b;
        T[i][k + 1] = -s * a + c * b;
    }
}

EigenResult eigendecompose_tridiagonal(double* const* A, int n, double tol) {
    EigenResult failed = {NULL, NULL, 0, 0};
    
    double** T = copy_matrix(A, n, n);
    if (T == NULL) return failed;
    
    double** V = identity(n);
    if (V == NULL) {
        free_matrix(T, n);
        return failed;
    }

    int max_iter_per_value = 1000;
    double* ks_c = malloc((n - 1) * sizeof(double));
    double* ks_s = malloc((n - 1) * sizeof(double));
    int*    ks_k = malloc((n - 1) * sizeof(int));
    
    if (ks_c == NULL || ks_s == NULL || ks_k == NULL) {
        free(ks_c);
        free(ks_s);
        free(ks_k);
        free_matrix(T, n);
        free_matrix(V, n);
        return failed;
    }

    int m = n;
    while (m > 1) {
        for (int iter = 0; iter < max_iter_per_value; iter++) {
            double mu = wilkinson_shift_tridiagonal(T[m-2][m-2], T[m-2][m-1], T[m-1][m-1]);

            for (int i = 0; i < m; i++) T[i][i] -= mu;

            int count = 0;
            for (int k = 0; k < m - 1; k++) {
                double a = T[k][k], b = T[k + 1][k];
                double r = hypot(a, b);
                double c = (r < 1e-300) ? 1.0 : a / r;
                double s = (r < 1e-300) ? 0.0 : b / r;
                ks_k[count] = k; ks_c[count] = c; ks_s[count] = s; count++;
                givens_left(T, k, c, s, m);
            }
            for (int idx = 0; idx < count; idx++) {
                givens_right(T, ks_k[idx], ks_c[idx], ks_s[idx], m);
                givens_right(V, ks_k[idx], ks_c[idx], ks_s[idx], n);
            }

            for (int i = 0; i < m; i++) T[i][i] += mu;

            if (fabs(T[m-1][m-2]) < tol) break;
        }
        m--;
    }

    free(ks_c); free(ks_s); free(ks_k);

    double* eigenvalues = diag(T, n);
    if (eigenvalues == NULL) {
        free_matrix(T, n);
        free_matrix(V, n);
        return failed;
    }
    
    free_matrix(T, n);

    double** eigenvectors_rows = transpose(V, n, n);
    if (eigenvectors_rows == NULL) {
        free(eigenvalues);
        free_matrix(V, n);
        return failed;
    }
    
    free_matrix(V, n);

    return (EigenResult){
        .eigenvalues = eigenvalues,
        .eigenvectors = eigenvectors_rows,
        .n = n,
        .count = n
    };
}

void free_eigen_result(EigenResult* r) {
    if (r->eigenvectors) {
        free_matrix(r->eigenvectors, r->count);
    }
    free(r->eigenvalues);
    *r = (EigenResult){0};
}

// The below commented-out function calculates eigenvectors using a null space,
// but was too unreliable due to eigenvalues being approximate. This is now rolled
// into the combined eigendecompose() function.
//
// EigenResult eigenvectors(double* const* A, int n, double tol) {
//     EigenResult failed = {0};
//     double* lambdas = eigenvalues(A, n, tol);
//     if (lambdas == NULL) return failed;

//     Node* eigenvector_list = create_linked_list();
//     if (eigenvector_list == NULL) {
//         free(lambdas);
//         return failed;
//     }

//     for (int i = 0; i < n; i++) {
//         double** in = identity(n);
//         double** lambdaIn = in ? multiply_matrix_by_scalar(in, n, n, lambdas[i]) : NULL;
//         double** system = lambdaIn ? subtract_matrices(lambdaIn, n, n, A, n, n) : NULL;
//         if (lambdaIn) free_matrix(lambdaIn, n);

//         double** transposed = system ? transpose(system, n, n) : NULL;
//         if (system) free_matrix(system, n);

//         double** AIm = transposed ? partition(n, transposed, n, in, n) : NULL;
//         if (transposed) free_matrix(transposed, n);
//         if (in) free_matrix(in, n);

//         if (AIm == NULL) {
//             free(lambdas);
//             free_linked_list(eigenvector_list, free);
//             return failed;
//         }

//         EliminationResult result = reduce(AIm, n, n * 2, n, sqrt(tol));
//         free_matrix(AIm, n);

//         int rank = length(result.pivots);
//         for (int j = 0; j < n; j++) {
//         }

//         bool alloc_failed = false;
//         for (int row = n - 1; row >= rank; row--) {
//             double* eigenvector = (double*)malloc(n * sizeof(double));
//             if (eigenvector == NULL) { alloc_failed = true; break; }
//             for (int j = n; j < n * 2; j++) {
//                 eigenvector[j - n] = result.rref[row][j];
//             }
//             push(eigenvector_list, eigenvector);
//         }

//         free_matrix(result.rref, n);
//         free_linked_list(result.pivots, free);

//         if (alloc_failed) {
//             free(lambdas);
//             free_linked_list(eigenvector_list, free);
//             return failed;
//         }
//     }

//     int num_eigenvectors = length(eigenvector_list);
//     double** result = (double**)malloc(num_eigenvectors * sizeof(double*));
//     if (result == NULL) {
//         free(lambdas);
//         free_linked_list(eigenvector_list, free);
//         return failed;
//     }

//     Node *current_eigenvector = eigenvector_list;
//     for (int i = 0; i < num_eigenvectors; i++) {
//         result[i] = (double*)current_eigenvector->data;
//         Node *next = current_eigenvector->next;
//         current_eigenvector = next;
//     }
//     free_linked_list(eigenvector_list, NULL);

//     return (EigenResult){
//         .eigenvalues = lambdas,
//         .eigenvectors = result,
//         .n = n,
//         .count = num_eigenvectors
//     };
// }
