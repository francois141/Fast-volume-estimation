#include <assert.h>
#include <stdio.h>

#include <armadillo>

#include "test_framework.h"

// extern "C" {
#include "linalg_optimized.h"
#include "polytope_optimized.h"
// }

/* MODIFY THIS IF MORE TESTS ARE ADDED */
#define TESTS_ARRAY_LEN 3

int perform_matrix_multiplication(int m, int n, int k) {
    int sizeof_a = m * k;
    int sizeof_b = k * n;
    int sizeof_c = m * n;

    double* A = (double*)malloc(sizeof(double) * sizeof_a);
    double* B = (double*)malloc(sizeof(double) * sizeof_b);
    double* C = (double*)malloc(sizeof(double) * sizeof_c);

    fill_with_random(A, sizeof_a);
    fill_with_random(B, sizeof_b);
    fill_with_zero(C, sizeof_c);

    arma::mat A_arma(m, k), B_arma(k, n);

    fill_with_same_values(A, &A_arma, m, k);
    fill_with_same_values(B, &B_arma, k, n);

    matrix_multiplication_wrapper(m, n, k, 1, 0, A, B, C);
    arma::mat C_arma = A_arma * B_arma;

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            // assert_same(C[i *n + j], C_arma(i,j), "Expected: %lf. Got: %lf\n", C[i * n + j], C_arma(i, j));
            assert_same(C[i *n + j], C_arma(i,j), "");
        }
    }

    free(A);
    free(B);
    free(C);
    return 0;
}

int test_matrix_multiplication_k_dominant() {
    int m = 4;
    int n = 1;
    int k = 6;
    return perform_matrix_multiplication(m, n, k);
}

int test_matrix_multiplication_n_dominant() {
    int m = 48;
    int n = 71;
    int k = 57;
    return perform_matrix_multiplication(m, n, k);
}

int test_matrix_multiplication_m_dominant() {
    int m = 443;
    int n = 243;
    int k = 105;
    return perform_matrix_multiplication(m, n, k);
}

int main(void) {
    /* Run all tests
     *  Current version simply calls the functions of the tests
     */

    srand((unsigned)time(NULL));

    int (*tests[TESTS_ARRAY_LEN])() = {&test_matrix_multiplication_m_dominant,
                                       &test_matrix_multiplication_n_dominant,
                                       &test_matrix_multiplication_k_dominant};
    for (int i = 0; i < TESTS_ARRAY_LEN; ++i) {// Make color red for any failure
        check_test(tests[i](),"Test %d succeeded\n", "Test %d failed (check messages)\n", i+1);
    }
    return 0;
}