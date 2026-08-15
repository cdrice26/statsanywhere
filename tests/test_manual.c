#define UNITY_INCLUDE_DOUBLE
#include "unity.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "operations.h"
#include "special_functions.h"

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

