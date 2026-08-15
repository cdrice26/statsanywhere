# API Documentation

## binomial.h

### `double binompmf(int n, int k, double p)`

**Binomial probability mass function.**

**Parameters:**
- `n` — Number of trials
- `k` — Number of successes
- `p` — Probability of success

**Returns:** Probability mass P(X = k)

### `double binomcdf(int n, double k, double p)`

**Binomial cumulative distribution function (lower tail).**

**Parameters:**
- `n` — Number of trials
- `k` — Threshold (may be non-integer; floor(k) is used)
- `p` — Probability of success

**Returns:** Cumulative probability P(X <= k)

## chi_squared.h

### `double X2pdf(double x, int df)`

**Chi-squared probability density function.**

**Parameters:**
- `x` — Point to evaluate (must be >= 0)
- `df` — Degrees of freedom

**Returns:** PDF value

### `double X2cdf(double x, int df)`

**Chi-squared cumulative distribution function.**

**Parameters:**
- `x` — Upper limit of integration
- `df` — Degrees of freedom

**Returns:** CDF value P(X <= x)

## f.h

### `double fpdf(double x, double d1, double d2)`

**F-distribution probability density function.**

**Parameters:**
- `x` — Point to evaluate
- `d1` — Numerator degrees of freedom
- `d2` — Denominator degrees of freedom

**Returns:** PDF value

### `double fcdf(double x, double d1, double d2)`

**F-distribution cumulative distribution function (lower tail).**

**Parameters:**
- `x` — Upper limit of integration
- `d1` — Numerator degrees of freedom
- `d2` — Denominator degrees of freedom

**Returns:** CDF value

## linked_list.h

### `Node* create_linked_list(void)`

**Creates an empty linked list head node.**

**Returns:** Pointer to newly-allocated head Node, or NULL on allocation failure

### `void free_linked_list(Node *head, void (*free_data)`

**Frees a linked list and optionally its element data.**

**Parameters:**
- `head` — Pointer to the head of the list (may be NULL)
- `free_data` — Optional callback to free each element's data; pass

### `int length(Node *head)`

**Returns the number of stored elements in the list.**

**Parameters:**
- `head` — Pointer to the head of the list

**Returns:** Number of elements stored in the list

## matrices.h

### `double** alloc_matrix(int m, int n)`

**Allocate an m-by-n matrix.**

**Parameters:**
- `m` — Number of rows
- `n` — Number of columns

**Returns:** Newly-allocated matrix (double**), or NULL on allocation failure

### `void free_matrix(double** A, int m)`

**Free an m-row matrix previously returned by alloc_matrix.**

**Parameters:**
- `A` — Matrix to free
- `m` — Number of rows in A

### `double** copy_matrix(double* const* A, int m, int n)`

**Make a deep copy of an m-by-n matrix.**

**Parameters:**
- `A` — Source matrix
- `m` — Number of rows
- `n` — Number of columns

**Returns:** Newly-allocated copy of A, or NULL on allocation failure

### `double** transpose(double* const* A, int m, int n)`

**Transpose an m-by-n matrix.**

**Parameters:**
- `A` — Source matrix
- `m` — Number of rows in A
- `n` — Number of columns in A

**Returns:** Transposed matrix (n-by-m), or NULL on allocation failure

### `double** gram_schmidt(double* const* A, int n)`

**Perform Gram–Schmidt orthonormalization on n row-vectors.**

**Parameters:**
- `A` — Array of n input row-vectors
- `n` — Dimension (number of rows and length of each row)

**Returns:** Newly-allocated orthonormal rows (Q), or NULL on allocation failure

### `QR_Decomposition QR_decompose(double* const* A, int n)`

**Compute the QR decomposition of an n-by-n matrix.**

**Parameters:**
- `A` — Input n-by-n matrix
- `n` — Matrix dimension

**Returns:** QR_Decomposition with .Q and .R set, or {NULL,NULL} on failure

### `double* diag(double* const* A, int n)`

**Extract the diagonal of an n-by-n matrix.**

**Parameters:**
- `A` — Input matrix
- `n` — Dimension

**Returns:** Newly-allocated array of diagonal elements, or NULL on failure

### `double** add_matrices(double* const* A, int ma, int na, double* const* B, int mb, int nb)`

**Add two matrices of identical shape.**

**Parameters:**
- `A` — First addend
- `ma` — Rows in A
- `na` — Columns in A
- `B` — Second addend
- `mb` — Rows in B
- `nb` — Columns in B

**Returns:** Newly-allocated sum matrix, or NULL on error

### `double** subtract_matrices(double* const* A, int ma, int na, double* const* B, int mb, int nb)`

**Subtract matrix B from A (A - B) for matching shapes.**

**Parameters:**
- `A` — Minuend matrix
- `ma` — Rows in A
- `na` — Columns in A
- `B` — Subtrahend matrix
- `mb` — Rows in B
- `nb` — Columns in B

**Returns:** Newly-allocated difference matrix, or NULL on error

### `double** multiply(double* const* A, int ma, int na, double* const* B, int mb, int nb)`

**Multiply two matrices A (ma-by-na) and B (mb-by-nb).**

**Parameters:**
- `A` — Left matrix
- `ma` — Rows in A
- `na` — Columns in A (must equal mb)
- `B` — Right matrix
- `mb` — Rows in B
- `nb` — Columns in B

**Returns:** Newly-allocated product matrix, or NULL on error

### `double** multiply_matrix_by_scalar(double* const* A, int rows, int cols, double s)`

**Multiply every element of a matrix by a scalar.**

**Parameters:**
- `A` — Source matrix
- `rows` — Number of rows
- `cols` — Number of columns
- `s` — Scalar multiplier

**Returns:** Newly-allocated scaled matrix, or NULL on allocation failure

### `double** identity(int n)`

**Create an n-by-n identity matrix.**

**Parameters:**
- `n` — Dimension

**Returns:** Newly-allocated identity matrix, or NULL on allocation failure

### `double** partition(int m, double* const* A, int an, double* const* B, int bn)`

**Horizontally concatenate two m-row matrices A and B.**

**Parameters:**
- `m` — Number of rows (must match for A and B)
- `A` — Left block matrix with an columns
- `an` — Number of columns in A
- `B` — Right block matrix with bn columns
- `bn` — Number of columns in B

**Returns:** Newly-allocated concatenated matrix, or NULL on allocation failure

### `EliminationResult reduce(double* const* A, int m, int n, int pivot_col_limit, double tol)`

**Compute row-reduced echelon form (RREF) with pivot tracking.**

**Parameters:**
- `A` — Input m-by-n matrix (not modified)
- `m` — Number of rows
- `n` — Number of columns
- `pivot_col_limit` — Maximum column index to consider for pivots
- `tol` — Tolerance for considering a value non-zero

**Returns:** EliminationResult with rref, dimensions, and pivot list; on failure rref will be NULL

### `EigenResult eigendecompose(double* const* A, int n, double tol)`

**Compute eigenvalues and eigenvectors (QR algorithm).**

**Parameters:**
- `A` — Input n-by-n matrix
- `n` — Dimension
- `tol` — Convergence tolerance for off-diagonal entries

**Returns:** EigenResult with eigenvalues and eigenvectors; on failure fields are NULL

### `EigenResult eigendecompose_tridiagonal(double* const* A, int n, double tol)`

**Compute eigen-decomposition for a symmetric tridiagonal matrix.**

**Parameters:**
- `A` — Input n-by-n tridiagonal matrix
- `n` — Dimension
- `tol` — Convergence tolerance

**Returns:** EigenResult with eigenvalues and eigenvectors; on failure fields are NULL

### `void free_eigen_result(EigenResult* r)`

**Free memory held by an EigenResult.**

**Parameters:**
- `r` — Pointer to EigenResult to free/clear

## normal.h

### `double normalpdf_std(double z)`

**Standard normal probability density function (mean 0, sd 1).**

**Parameters:**
- `z` — Standardized coordinate

**Returns:** PDF value

### `double normalcdf_std(double lower, double upper)`

**Standard normal probability between two z-values.**

**Parameters:**
- `lower` — Lower z-value
- `upper` — Upper z-value

**Returns:** Probability mass between lower and upper

### `double normalcdf(double lower, double upper, double mu, double sigma)`

**Normal CDF between two values for N(mu, sigma^2).**

**Parameters:**
- `lower` — Lower bound
- `upper` — Upper bound
- `mu` — Mean
- `sigma` — Standard deviation

**Returns:** Probability P(lower <= X <= upper)

## operations.h

### `double factorial(unsigned int n)`

**Compute the factorial of n.**

**Parameters:**
- `n` — Non-negative integer

**Returns:** n! as a double

⚠️ **Warning:** the implementation does not perform overflow checks for large n.

### `double combination(unsigned int n, unsigned int r)`

**Compute the binomial coefficient C(n, r).**

**Parameters:**
- `n` — Total number of items
- `r` — Number of selections

**Returns:** Binomial coefficient C(n, r) as a double

## poisson.h

### `double poissonpmf(int k, double lambda)`

**Poisson probability mass function.**

**Parameters:**
- `k` — Number of occurrences (integer)
- `lambda` — Rate parameter (>= 0)

**Returns:** Probability mass P(K = k)

### `double poissoncdf(double k, double lambda)`

**Poisson cumulative distribution function (lower tail).**

**Parameters:**
- `k` — Threshold (may be non-integer; floor(k) is used)
- `lambda` — Rate parameter

**Returns:** Cumulative probability P(K <= k)

## quadrature_weights.h

### `double quadrature_estimate(int n, double* a, double* b, double mu_0, double (*f)`

**Estimate an integral using Gaussian quadrature from Jacobi matrix.**

**Parameters:**
- `n` — Number of quadrature nodes
- `a` — Array of length n containing diagonal generators
- `b` — Array of length n-1 containing squared off-diagonal generators
- `mu_0` — Zeroth moment / scaling factor for weight computation
- `f` — Integrand function of a node and user context
- `ctx` — User-provided context forwarded to f

**Returns:** Estimated integral value, or negative on failure

## special_functions.h

### `double gamma(double z)`

**Compute the Gamma function via quadrature.**

**Parameters:**
- `z` — Argument to the Gamma function

**Returns:** Approximated Gamma(z)

### `double lower_incomplete_gamma(double s, double x)`

**Compute the lower incomplete gamma function.**

**Parameters:**
- `s` — Shape parameter
- `x` — Integration upper limit parameter

**Returns:** Approximated lower incomplete gamma value

### `double incomplete_beta(double x, double a, double b)`

**Compute the incomplete beta function.**

**Parameters:**
- `x` — Upper integration limit parameter
- `a` — First shape parameter
- `b` — Second shape parameter

**Returns:** Approximated incomplete beta value (unnormalized)

### `double regularized_incomplete_beta(double x, double a, double b)`

**Compute the regularized incomplete beta function.**

**Parameters:**
- `x` — Upper integration limit parameter
- `a` — First shape parameter
- `b` — Second shape parameter

**Returns:** Regularized incomplete beta in [0,1]

### `double beta(double z1, double z2)`

**Compute the Beta function via Gamma functions.**

**Parameters:**
- `z1` — First argument
- `z2` — Second argument

**Returns:** Beta(z1, z2)

### `double erf(double z)`

**Compute the error function erf(z).**

**Parameters:**
- `z` — Argument to erf

**Returns:** Approximated erf(z)

## t.h

### `double tpdf(double t, int df)`

**Student's t-distribution probability density function.**

**Parameters:**
- `t` — Point to evaluate
- `df` — Degrees of freedom

**Returns:** PDF value

### `double tcdf(double t, int df)`

**Student's t-distribution cumulative distribution function.**

**Parameters:**
- `t` — Point to evaluate
- `df` — Degrees of freedom

**Returns:** CDF value

⚠️ **Warning:** Not quite as precise as some other functions due to multiple quadratures being used

## vectors.h

### `double* add(const double* v1, const double* v2, int n)`

**Adds two vectors and returns the resultant.**

**Parameters:**
- `v1` — The first vector
- `v2` — The second vector

**Returns:** A newly-allocated resultant vector.

### `double* subtract(const double* v1, const double* v2, int n)`

**Subtracts one vector from another.**

**Parameters:**
- `v1` — Minuend vector
- `v2` — Subtrahend vector
- `n` — Number of elements in each vector

**Returns:** Newly-allocated result vector, or NULL on allocation failure

### `double dot_product(const double* v1, const double* v2, int n)`

**Computes the dot product of two vectors.**

**Parameters:**
- `v1` — First input vector
- `v2` — Second input vector
- `n` — Number of elements in each vector

**Returns:** Scalar dot product

### `double* multiply_scalar(const double* v, int n, double k)`

**Multiplies a vector by a scalar.**

**Parameters:**
- `v` — Input vector
- `n` — Number of elements in v
- `k` — Scalar multiplier

**Returns:** Newly-allocated scaled vector, or NULL on allocation failure

### `double* projection(const double* of, const double* onto, int n)`

**Projects one vector onto another.**

**Parameters:**
- `of` — Vector being projected
- `onto` — Vector to project onto
- `n` — Number of elements in each vector

**Returns:** Newly-allocated projection vector, or NULL on allocation failure

### `double magnitude(const double* v, int n)`

**Computes the Euclidean magnitude (length) of a vector.**

**Parameters:**
- `v` — Input vector
- `n` — Number of elements in v

**Returns:** Euclidean magnitude

### `double* copy_vector(const double* v, int n)`

**Creates a copy of a vector.**

**Parameters:**
- `v` — Input vector
- `n` — Number of elements in v

**Returns:** Newly-allocated copy of v, or NULL on allocation failure

