#include "vectors.h"

/**
 * @brief Adds two vectors and returns the resultant.
 *
 * Returns the resultant of adding the two vectors. Caller is responsible for freeing resultant.
 *
 * @param v1 The first vector
 * @param v2 The second vector
 *
 * @return A newly-allocated resultant vector.
 */
double* add(const double* v1, const double* v2, int n) {
    double* result = (double*)malloc(n * sizeof(double));
    if (result == NULL) return NULL;
    for (int i = 0; i < n; i++) {
        result[i] = v1[i] + v2[i];
    }
    return result;
}

/**
 * @brief Negates a vector (element-wise).
 *
 * Allocates and returns a new vector where each element is the negation of
 * the corresponding element in v. Caller is responsible for freeing the
 * returned array.
 *
 * @param v Input vector
 * @param n Number of elements in v
 * @return Newly-allocated negated vector, or NULL on allocation failure
 */
double* negate(const double* v, int n) {
    double* result = (double*)malloc(n * sizeof(double));
    if (result == NULL) return NULL;
    for (int i = 0; i < n; i++) {
        result[i] = -v[i];
    }
    return result;
}

/**
 * @brief Subtracts one vector from another.
 *
 * Computes v1 - v2 and returns a newly-allocated vector with the result.
 * This function uses negate() and add() internally. Caller must free the
 * returned array.
 *
 * @param v1 Minuend vector
 * @param v2 Subtrahend vector
 * @param n Number of elements in each vector
 * @return Newly-allocated result vector, or NULL on allocation failure
 */
double* subtract(const double* v1, const double* v2, int n) {
    double* nv2 = negate(v2, n);
    if (nv2 == NULL) return NULL;
    double* result = add(v1, nv2, n);
    if (result == NULL) return NULL;
    free(nv2);
    return result;
}

/**
 * @brief Computes the dot product of two vectors.
 *
 * Returns the scalar dot product sum(v1[i] * v2[i]) for i in [0, n).
 *
 * @param v1 First input vector
 * @param v2 Second input vector
 * @param n Number of elements in each vector
 * @return Scalar dot product
 */
double dot_product(const double* v1, const double* v2, int n) {
    double result = 0;
    for (int i = 0; i < n; i++) {
        result += v1[i] * v2[i];
    }
    return result;
}

/**
 * @brief Multiplies a vector by a scalar.
 *
 * Allocates and returns a new vector where each element is v[i] * k.
 * Caller is responsible for freeing the returned array.
 *
 * @param v Input vector
 * @param n Number of elements in v
 * @param k Scalar multiplier
 * @return Newly-allocated scaled vector, or NULL on allocation failure
 */
double* multiply_scalar(const double* v, int n, double k) {
    double* result = (double*)malloc(n * sizeof(double));
    if (result == NULL) return NULL;
    for (int i = 0; i < n; i++) {
        result[i] = v[i] * k;
    }
    return result;
}

/**
 * @brief Projects one vector onto another.
 *
 * Returns a newly-allocated vector representing the projection of 'of' onto
 * 'onto', computed as onto * (dot(onto, of) / dot(onto, onto)). Caller must
 * free the returned array.
 *
 * @param of Vector being projected
 * @param onto Vector to project onto
 * @param n Number of elements in each vector
 * @return Newly-allocated projection vector, or NULL on allocation failure
 */
double* projection(const double* of, const double* onto, int n) {
    return multiply_scalar(onto, n, dot_product(onto, of, n) / dot_product(onto, onto, n));
}

/**
 * @brief Computes the Euclidean magnitude (length) of a vector.
 *
 * Returns sqrt(dot(v, v)).
 *
 * @param v Input vector
 * @param n Number of elements in v
 * @return Euclidean magnitude
 */
double magnitude(const double* v, int n) {
    return sqrt(dot_product(v, v, n));
}

/**
 * @brief Creates a copy of a vector.
 *
 * Allocates and returns a new vector containing the same elements as v.
 * Caller is responsible for freeing the returned array.
 *
 * @param v Input vector
 * @param n Number of elements in v
 * @return Newly-allocated copy of v, or NULL on allocation failure
 */
double* copy_vector(const double* v, int n) {
    return multiply_scalar(v, n, 1.0);
}

/**
 * @brief Finds the maximum element in a vector.
 *
 * Scans the vector and returns the largest value.
 *
 * @param v Input vector
 * @param n Number of elements in v (must be >= 1)
 * @return Maximum element value
 */
double max(const double* v, int n) {
    double max = v[0];
    for (int i = 1; i < n; i++) {
        if (v[i] > max) max = v[i];
    }
    return max;
}
