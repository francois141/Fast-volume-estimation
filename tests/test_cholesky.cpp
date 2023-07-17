//
// Created by francois on 23.04.23.
//

#include <stdio.h>
#include <armadillo>
#include "test_framework.h"

// extern "C" {
#include "linalg_optimized.h"
#include "polytope_optimized.h"
// }


int smoke_test_cholesky_armadillo() {

    arma::mat A_arma(2, 2);

    A_arma(0,0) = 1;
    A_arma(0,1) = 1;
    A_arma(1,0) = 1;
    A_arma(1,1) = 11;

    arma::mat R = arma::chol(A_arma);

    double expect[2][2] = {{1,1},{0,3.162278}};

    for(int i = 0; i < 2;i++) {
        for(int j = 0; j < 2;j++) {
            int idx = 2*i + j;
            // assert_same(expect[i][j], R(i,j),"Expected result: %lf. Result armadillo: %lf\n", expect[i][j], R(i,j));
            assert_same(expect[i][j], R(i,j),"");
        }
    }

    return 0;
}

int test_wrapper() {

    #define MAX_SIZE 16
    struct test_entry {
        int n;
        double arr[MAX_SIZE];
        double expect[MAX_SIZE];
    };

    const unsigned int nb_tests = 4;

    struct test_entry tests[nb_tests];

    tests[0] = (struct test_entry){
            .n = 4,
            .arr = {1,1,1,1,1,5,5,5,1,5,14,14,1,5,14,15},
            .expect = {1,1,1,1,0,2,2,2,0,0,3,3,0,0,0,1},
    };

    tests[1] = (struct test_entry){
            .n = 2,
            .arr = {4,3,3,4},
            .expect = {2,1.5,0,1.323},
    };

    tests[2] = (struct test_entry){
            .n = 2,
            .arr = {1,1,1,11},
            .expect = {1,1,0,3.162},
    };

    tests[3] = (struct test_entry){
            .n = 3,
            .arr = {2,1,1,1,2,1,1,1,2},
            .expect = {1.414214,0.707107,0.707107,0,1.224745,0.408248, 0,0,1.154701},
    };

    for(int test_idx = 0; test_idx < nb_tests;test_idx++) {

        const unsigned int matrix_size = tests[test_idx].n;

        arma::mat A_armadillo(matrix_size, matrix_size);
        fill_with_same_values(tests[test_idx].arr, &A_armadillo, matrix_size, matrix_size);

        int error_code = cholesky_decomposition(matrix_size, tests[test_idx].arr);
        // assert_zero(error_code, "Error code should be equal 0 : %d\n", error_code);
        assert_zero(error_code, "");

        for(int i = 0; i < matrix_size;i++) {
            for(int j = 0; j < matrix_size;j++) {
                int idx = matrix_size*i + j;
                assert_same(tests[test_idx].arr[idx],tests[test_idx].expect[idx], "");
            }
            // printf("\n");
        }

        arma::mat R = arma::chol(A_armadillo);

        for(int i = 0; i < matrix_size;i++) {
            for(int j = 0; j < matrix_size;j++) {
                int idx = matrix_size*i + j;
                // assert_same(tests[test_idx].arr[idx], R(i,j),"Result wrapper: %lf. Result armadillo: %lf\n", tests[test_idx].arr[idx], R(i,j));
                assert_same(tests[test_idx].arr[idx], R(i,j),"");
            }
        }
    }

    return 0;
}

int main(void) {
    srand((unsigned)time(NULL));

    const unsigned int nb_tests = 2;

    int (*tests[nb_tests])() = {&test_wrapper, &smoke_test_cholesky_armadillo};

    for (int i = 0; i < nb_tests; ++i) {// Make color red for any failure
        check_test(tests[i](),"Test %d succeeded\n", "Test %d failed (check messages)\n", i+1);
    }

    return 0;
}