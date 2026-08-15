#include "unity.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "operations.h"
#include "special_functions.h"
#include "matrices.h"
#include "linked_list.h"

static const double EPS = 1e-6;

static int get_expected(const char *name, const char *args, double *out) {
        char line[256];
        if (!name || !out) return 0;
        if (!args) args = "";

#ifdef REFERENCE_FILE_PATH
        const char *ref_path = REFERENCE_FILE_PATH;
#else
        const char *ref_path = "tests/reference_values.txt";
#endif
        FILE *ref = fopen(ref_path, "r");
        if (!ref) {
            char msg[256];
            snprintf(msg, sizeof(msg), "Missing %s — run: cmake --build build --target generate_reference && cmake --build build --target test_runner", ref_path);
            TEST_FAIL_MESSAGE(msg);
            return 0;
        }

        while (fgets(line, sizeof(line), ref)) {
            char *last_space = strrchr(line, ' ');
            if (!last_space) continue;
            double val = atof(last_space + 1);
            *last_space = '\0';
            if (strcmp(line, name) == 0 || strstr(line, name) == line) {
                if (args[0] == '\0' || strstr(line, args) != NULL) {
                    *out = val;
                    fclose(ref);
                    return 1;
                }
            }
        }
        fclose(ref);
        TEST_FAIL_MESSAGE("Reference value not found in tests/reference_values.txt");
        return 0;
}

/* Helper: compute the maximum absolute difference between two matrices */
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

/* setUp/tearDown are provided by the test runner */

void test_factorial_and_combination(void) {
    double expected;
    if (!get_expected("factorial", "9", &expected)) { TEST_IGNORE_MESSAGE("reference generator unavailable"); }
    TEST_ASSERT_EQUAL_DOUBLE(expected, factorial(9));

    if (!get_expected("combination", "6 4", &expected)) { TEST_IGNORE_MESSAGE("reference generator unavailable"); }
    TEST_ASSERT_EQUAL_DOUBLE(expected, combination(6, 4));
}

void test_gamma_and_incomplete(void) {
    double expected;
    if (!get_expected("gamma", "5.0", &expected)) { TEST_IGNORE_MESSAGE("reference generator unavailable"); }
    TEST_ASSERT_DOUBLE_WITHIN(EPS, expected, gamma(5.0));

    if (!get_expected("gamma", "6.7", &expected)) { TEST_IGNORE_MESSAGE("reference generator unavailable"); }
    TEST_ASSERT_DOUBLE_WITHIN(1e-5, expected, gamma(6.7));

    if (!get_expected("lower_incomplete_gamma", "1 5", &expected)) { TEST_IGNORE_MESSAGE("reference generator unavailable"); }
    TEST_ASSERT_DOUBLE_WITHIN(1e-4, expected, lower_incomplete_gamma(1, 5));

    if (!get_expected("lower_incomplete_gamma", "3 5.5", &expected)) { TEST_IGNORE_MESSAGE("reference generator unavailable"); }
    TEST_ASSERT_DOUBLE_WITHIN(1e-4, expected, lower_incomplete_gamma(3, 5.5));
}

void test_beta_and_incomplete_beta(void) {
    double expected;
    if (!get_expected("beta", "3 2", &expected)) { TEST_IGNORE_MESSAGE("reference generator unavailable"); }
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, expected, beta(3, 2));

    if (!get_expected("beta", "4 10", &expected)) { TEST_IGNORE_MESSAGE("reference generator unavailable"); }
    TEST_ASSERT_DOUBLE_WITHIN(5e-7, expected, beta(4, 10));

    if (!get_expected("incomplete_beta", "0.4 1 3", &expected)) { TEST_IGNORE_MESSAGE("reference generator unavailable"); }
    TEST_ASSERT_DOUBLE_WITHIN(1e-4, expected, incomplete_beta(0.4, 1, 3));

    if (!get_expected("incomplete_beta", "1 2 3", &expected)) { TEST_IGNORE_MESSAGE("reference generator unavailable"); }
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, expected, incomplete_beta(1, 2, 3));

    if (!get_expected("regularized_incomplete_beta", "0.4 1 3", &expected)) { TEST_IGNORE_MESSAGE("reference generator unavailable"); }
    TEST_ASSERT_DOUBLE_WITHIN(1e-4, expected, regularized_incomplete_beta(0.4, 1, 3));

    if (!get_expected("regularized_incomplete_beta", "1 2 3", &expected)) { TEST_IGNORE_MESSAGE("reference generator unavailable"); }
    TEST_ASSERT_DOUBLE_WITHIN(1e-8, expected, regularized_incomplete_beta(1, 2, 3));
}

void test_erf_values(void) {
    double expected;
    if (!get_expected("erf", "1", &expected)) { TEST_IGNORE_MESSAGE("reference generator unavailable"); }
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, expected, erf(1));

    if (!get_expected("erf", "0.5", &expected)) { TEST_IGNORE_MESSAGE("reference generator unavailable"); }
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, expected, erf(0.5));
}

void test_matrix_functions(void) {
    int n = 3;

    double** A = alloc_matrix(n, n);
    A[0][0] = 2;  A[0][1] = -1; A[0][2] = 4;
    A[1][0] = -1;   A[1][1] = 3; A[1][2] = 0;
    A[2][0] = 4;  A[2][1] = 0;  A[2][2] = 5;

    EliminationResult res = reduce(A, n, n, n, 1e-6);
    double** rref = res.rref;
    void** pivots = get_data(res.pivots);
    for (int i = 0; i < 3; i++) {
        int val = *(int*)pivots[i];
        TEST_ASSERT_EQUAL_INT(i, val);
    }
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 1, rref[0][0]);
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 1, rref[1][1]);
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 1, rref[2][2]);
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 0, rref[0][1]);
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 0, rref[1][2]);
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 0, rref[2][0]);
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 0, rref[0][2]);
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 0, rref[1][0]);
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 0, rref[2][1]);
    TEST_ASSERT_EQUAL_INT(3, res.m);
    TEST_ASSERT_EQUAL_INT(3, res.n);

    free_linked_list(res.pivots, NULL);
    free(pivots);
    free_matrix(res.rref, n);

    QR_Decomposition qr = QR_decompose(A, n);
    double** QR = multiply(qr.Q, n, n, qr.R, n, n);
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, 0, max_abs_diff(QR, A, n));
    free_matrix(qr.Q, n);
    free_matrix(qr.R, n);
    free_matrix(QR, n);
    free_matrix(A, n);
}
