#include "matrices.h"
#include "linked_list.h"

/**
 * @brief Allocate an m-by-n matrix.
 *
 * Allocates an array of m pointers to rows and m row arrays of length n.
 * Caller is responsible for freeing with free_matrix.
 *
 * @param m Number of rows
 * @param n Number of columns
 * @return Newly-allocated matrix (double**), or NULL on allocation failure
 */
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

/**
 * @brief Free an m-row matrix previously returned by alloc_matrix.
 *
 * Frees each row and then the row-pointer array.
 *
 * @param A Matrix to free
 * @param m Number of rows in A
 */
void free_matrix(double** A, int m) {
    for (int i = 0; i < m; i++) {
        free(A[i]);
    }
    free(A);
}

/**
 * @brief Make a deep copy of an m-by-n matrix.
 *
 * Allocates and returns a new matrix with the same contents as A. Caller
 * must free the returned matrix with free_matrix.
 *
 * @param A Source matrix
 * @param m Number of rows
 * @param n Number of columns
 * @return Newly-allocated copy of A, or NULL on allocation failure
 */
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

/**
 * @brief Transpose an m-by-n matrix.
 *
 * Allocates and returns the n-by-m transpose of A. Caller must free the
 * returned matrix with free_matrix.
 *
 * @param A Source matrix
 * @param m Number of rows in A
 * @param n Number of columns in A
 * @return Transposed matrix (n-by-m), or NULL on allocation failure
 */
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

/**
 * @brief Perform Gram–Schmidt orthonormalization on n row-vectors.
 *
 * Treats A as an array of n row vectors of length n and returns an array of
 * n orthonormal row vectors (Q). Caller must free the returned array and
 * its rows with free_matrix.
 *
 * @param A Array of n input row-vectors
 * @param n Dimension (number of rows and length of each row)
 * @return Newly-allocated orthonormal rows (Q), or NULL on allocation failure
 */
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

/**
 * @brief Compute the QR decomposition of an n-by-n matrix.
 *
 * Returns Q and R such that A = Q * R where Q is orthonormal and R is
 * upper-triangular. The returned QR_Decomposition contains allocated Q and R
 * which the caller must free with free_matrix.
 *
 * @param A Input n-by-n matrix
 * @param n Matrix dimension
 * @return QR_Decomposition with .Q and .R set, or {NULL,NULL} on failure
 */
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

/**
 * @brief Extract the diagonal of an n-by-n matrix.
 *
 * Allocates and returns an array of length n containing A[i][i]. Caller must
 * free the returned array.
 *
 * @param A Input matrix
 * @param n Dimension
 * @return Newly-allocated array of diagonal elements, or NULL on failure
 */
double* diag(double* const* A, int n) {
    double* result = (double*)malloc(n * sizeof(double));
    if (result == NULL) return NULL;
    for (int i = 0; i < n; i++) {
        result[i] = A[i][i];
    }
    return result;
}

/**
 * @brief Add two matrices of identical shape.
 *
 * Allocates and returns a new ma-by-na matrix equal to A + B. Returns NULL if
 * shapes mismatch or allocation fails. Caller must free the result.
 *
 * @param A First addend
 * @param ma Rows in A
 * @param na Columns in A
 * @param B Second addend
 * @param mb Rows in B
 * @param nb Columns in B
 * @return Newly-allocated sum matrix, or NULL on error
 */
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

/**
 * @brief Subtract matrix B from A (A - B) for matching shapes.
 *
 * Allocates and returns a new matrix with the element-wise difference. Caller
 * must free the returned matrix with free_matrix.
 *
 * @param A Minuend matrix
 * @param ma Rows in A
 * @param na Columns in A
 * @param B Subtrahend matrix
 * @param mb Rows in B
 * @param nb Columns in B
 * @return Newly-allocated difference matrix, or NULL on error
 */
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


/**
 * @brief Multiply two matrices A (ma-by-na) and B (mb-by-nb).
 *
 * Returns a newly-allocated ma-by-nb matrix equal to A * B. Returns NULL if
 * inner dimensions mismatch or allocation fails. Caller must free result.
 *
 * @param A Left matrix
 * @param ma Rows in A
 * @param na Columns in A (must equal mb)
 * @param B Right matrix
 * @param mb Rows in B
 * @param nb Columns in B
 * @return Newly-allocated product matrix, or NULL on error
 */
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

/**
 * @brief Multiply every element of a matrix by a scalar.
 *
 * Allocates and returns a new matrix where each element is A[i][j] * s.
 * Caller must free the returned matrix.
 *
 * @param A Source matrix
 * @param rows Number of rows
 * @param cols Number of columns
 * @param s Scalar multiplier
 * @return Newly-allocated scaled matrix, or NULL on allocation failure
 */
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

/**
 * @brief Create an n-by-n identity matrix.
 *
 * Allocates and returns an identity matrix with 1.0 on the diagonal and
 * 0.0 elsewhere. Caller must free the returned matrix.
 *
 * @param n Dimension
 * @return Newly-allocated identity matrix, or NULL on allocation failure
 */
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

/**
 * @brief Horizontally concatenate two m-row matrices A and B.
 *
 * Allocates and returns an m-by-(an+bn) matrix whose left block is A and
 * right block is B. Caller must free the returned matrix.
 *
 * @param m Number of rows (must match for A and B)
 * @param A Left block matrix with an columns
 * @param an Number of columns in A
 * @param B Right block matrix with bn columns
 * @param bn Number of columns in B
 * @return Newly-allocated concatenated matrix, or NULL on allocation failure
 */
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

/**
 * @brief Swap two rows of a matrix in-place.
 *
 * Exchanges the row pointers A[i] and A[j]. If i == j, nothing happens.
 *
 * @param A Matrix whose rows are swapped
 * @param i Index of first row
 * @param j Index of second row
 */
void swap_rows(double **A, int i, int j) {
    if (i == j) return;
    double *tmp = A[i];
    A[i] = A[j];
    A[j] = tmp;
}

/**
 * @brief Scale a row by a constant factor.
 *
 * Multiplies each element of row i by c across the given number of columns.
 *
 * @param A Matrix
 * @param cols Number of columns in A
 * @param i Row index to scale
 * @param c Scaling factor
 */
void scale_row(double **A, int cols, int i, double c) {
    for (int k = 0; k < cols; ++k) A[i][k] *= c;
}

/**
 * @brief Add a multiple of row j to row i (row_i += c * row_j).
 *
 * If c is zero the function returns early. Operates across the given number
 * of columns.
 *
 * @param A Matrix
 * @param cols Number of columns in A
 * @param i Destination row index
 * @param j Source row index
 * @param c Multiplier applied to row j before adding
 */
void add_multiple_row(double **A, int cols, int i, int j, double c) {
    if (c == 0.0) return;
    for (int k = 0; k < cols; ++k) A[i][k] += c * A[j][k];
}

/**
 * @brief Perform the back-substitution step of Gauss elimination.
 *
 * Given an EliminationResult containing rref and pivot positions, produces
 * a matrix with zeros above pivots by eliminating upward. Returns a newly
 * allocated matrix which the caller must free.
 *
 * @param result Pointer to EliminationResult containing rref, dimensions and pivots
 * @return Newly-allocated post-back-substitution matrix, or NULL on error
 */
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

/**
 * @brief Compute row-reduced echelon form (RREF) with pivot tracking.
 *
 * Performs forward elimination with partial pivoting up to pivot_col_limit
 * and returns an EliminationResult containing the rref matrix, its
 * dimensions, and a linked-list of pivot Points. The returned rref matrix is
 * allocated and must be freed by the caller. The pivot list must also be
 * freed with free_linked_list.
 *
 * @param A Input m-by-n matrix (not modified)
 * @param m Number of rows
 * @param n Number of columns
 * @param pivot_col_limit Maximum column index to consider for pivots
 * @param tol Tolerance for considering a value non-zero
 * @return EliminationResult with rref, dimensions, and pivot list; on
 *         failure rref will be NULL
 */
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

/**
 * @brief Compute Wilkinson shift for the bottom-right 2x2 block.
 *
 * Returns the Wilkinson shift value used to accelerate QR iteration for
 * eigenvalue computation on the trailing m-by-m block.
 *
 * @param A Matrix containing the trailing block
 * @param n Dimension of the (square) matrix
 * @return Wilkinson shift scalar
 */
double wilkinson_shift(double** A, int n) {
    double a = A[n-2][n-2];
    double b = A[n-2][n-1];
    double c = A[n-1][n-1];
    double delta = (a - c) / 2.0;
    double sign = (delta >= 0) ? 1.0 : -1.0;
    if (delta == 0.0) sign = 1.0;
    return c - sign * (b * b) / (fabs(delta) + sqrt(delta * delta + b * b));
}

/**
 * @brief Compute eigenvalues and eigenvectors (QR algorithm).
 *
 * Uses shifted QR iterations to compute eigenvalues and eigenvectors of an
 * n-by-n matrix A. Returns EigenResult containing allocated eigenvalues
 * (array of length n) and eigenvectors (n row-vectors). Caller must free
 * these using free_eigen_result.
 *
 * @param A Input n-by-n matrix
 * @param n Dimension
 * @param tol Convergence tolerance for off-diagonal entries
 * @return EigenResult with eigenvalues and eigenvectors; on failure fields are NULL
 */
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

/**
 * @brief Wilkinson shift for a 2x2 tridiagonal block.
 *
 * Variant of the Wilkinson shift specialized for tridiagonal matrices where
 * the 2x2 block is [a b; b c]. This function is static/internal.
 *
 * @param a Top-left element
 * @param b Off-diagonal element
 * @param c Bottom-right element
 * @return Shift value
 */
static double wilkinson_shift_tridiagonal(double a, double b, double c) {
    double delta = (a - c) / 2.0;
    double sign = (delta >= 0.0) ? 1.0 : -1.0;
    if (delta == 0.0) sign = 1.0;
    double denom = fabs(delta) + sqrt(delta * delta + b * b);
    return (denom != 0.0) ? c - sign * (b * b) / denom : c;
}

/**
 * @brief Apply Givens rotation to rows k and k+1 (left multiplication).
 *
 * Rotates rows k and k+1 in-place over the first `width` columns using
 * cosine c and sine s. Intended for applying Givens rotations from the
 * left.
 *
 * @param T Matrix to modify
 * @param k Row index of the first row
 * @param c Cosine component of rotation
 * @param s Sine component of rotation
 * @param width Number of columns to rotate
 */
static void givens_left(double** T, int k, double c, double s, int width) {
    for (int j = 0; j < width; j++) {
        double a = T[k][j], b = T[k + 1][j];
        T[k][j]     =  c * a + s * b;
        T[k + 1][j] = -s * a + c * b;
    }
}

/**
 * @brief Apply Givens rotation to columns k and k+1 (right multiplication).
 *
 * Rotates columns k and k+1 in-place over the first `rows` rows using
 * cosine c and sine s. Intended for applying Givens rotations from the
 * right.
 *
 * @param T Matrix to modify
 * @param k Column index of the first column
 * @param c Cosine component of rotation
 * @param s Sine component of rotation
 * @param rows Number of rows to rotate
 */
static void givens_right(double** T, int k, double c, double s, int rows) {
    for (int i = 0; i < rows; i++) {
        double a = T[i][k], b = T[i][k + 1];
        T[i][k]     =  c * a + s * b;
        T[i][k + 1] = -s * a + c * b;
    }
}

/**
 * @brief Compute eigen-decomposition for a symmetric tridiagonal matrix.
 *
 * Uses Givens rotations and shifted QR tailored to tridiagonal matrices to
 * compute eigenvalues and eigenvectors. Returns allocated eigenvalues and
 * eigenvectors (rows) which must be freed by the caller via
 * free_eigen_result.
 *
 * @param A Input n-by-n tridiagonal matrix
 * @param n Dimension
 * @param tol Convergence tolerance
 * @return EigenResult with eigenvalues and eigenvectors; on failure fields are NULL
 */
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

/**
 * @brief Free memory held by an EigenResult.
 *
 * Frees the eigenvectors matrix and eigenvalues array inside r and resets
 * the structure to zero.
 *
 * @param r Pointer to EigenResult to free/clear
 */
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
