#define UNITY_INCLUDE_DOUBLE
#include "unity.h"

/* declare test functions from other translation units */
void test_binomial_examples(void);
void test_chi_squared_examples(void);
void test_f_distribution_examples(void);
void test_normal_examples(void);
void test_poisson_examples(void);
void test_student_t_examples(void);

void test_factorial_and_combination(void);
void test_gamma_and_incomplete(void);
void test_beta_and_incomplete_beta(void);
void test_erf_values(void);

void setUp(void) {}
void tearDown(void) {}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_binomial_examples);
    RUN_TEST(test_chi_squared_examples);
    RUN_TEST(test_f_distribution_examples);
    RUN_TEST(test_normal_examples);
    RUN_TEST(test_poisson_examples);
    RUN_TEST(test_student_t_examples);

    RUN_TEST(test_factorial_and_combination);
    RUN_TEST(test_gamma_and_incomplete);
    RUN_TEST(test_beta_and_incomplete_beta);
    RUN_TEST(test_erf_values);

    return UNITY_END();
}
