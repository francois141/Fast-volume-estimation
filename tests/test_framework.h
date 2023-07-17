#ifndef TEST_FRAMEWORK_H_
#define TEST_FRAMEWORK_H_
#include <iostream>
#include <vector>

#include "math.h"

using namespace std;

#define RESET "\033[0m"
#define BLACK "\033[30m" /* Black */
#define RED "\033[31m"   /* Red */
#define GREEN "\033[32m" /* Green */

#define printf_green(x, ...) printf(GREEN x RESET, ##__VA_ARGS__)
#define printf_red(x, ...) printf(RED x RESET, ##__VA_ARGS__)

#define EPS 1e-3

#define assert_same(x, y, msg, ...)                           \
    if (((isinf(x) || isnan(x)) && (isinf(y) || isnan(y))) || \
        check_relative_difference(x, y)) {                    \
        printf_green(msg, ##__VA_ARGS__);                     \
    } else {                                                  \
        printf_red(msg, ##__VA_ARGS__);                       \
        exit(1);                                              \
    }

#define assert_not_same(x, y, msg, ...)     \
    if (!check_relative_difference(x, y)) { \
        printf_green(msg, ##__VA_ARGS__);   \
    } else {                                \
        printf_red(msg, ##__VA_ARGS__);     \
        exit(1);                            \
    }

#define assert_not_null(x, msg, ...)      \
    if (x != NULL) {                      \
        printf_green(msg, ##__VA_ARGS__); \
    } else {                              \
        printf_red(msg, ##__VA_ARGS__);   \
        exit(1);                          \
    }

#define assert_null(x, msg, ...)          \
    if (x == NULL) {                      \
        printf_green(msg, ##__VA_ARGS__); \
    } else {                              \
        printf_red(msg, ##__VA_ARGS__);   \
        exit(1);                          \
    }

#define assert_not_zero(x, msg, ...)      \
    if (x != 0) {                         \
        printf_green(msg, ##__VA_ARGS__); \
    } else {                              \
        printf_red(msg, ##__VA_ARGS__);   \
        exit(1);                          \
    }
#define assert_zero(x, msg, ...)          \
    if (x == 0) {                         \
        printf_green(msg, ##__VA_ARGS__); \
    } else {                              \
        printf_red(msg, ##__VA_ARGS__);   \
        exit(1);                          \
    }

#define check_test(ret_code, success_msg, fail_msg, ...) \
    if (ret_code == 0) {                                 \
        printf_green(success_msg, ##__VA_ARGS__);        \
    } else {                                             \
        printf_red(fail_msg, ##__VA_ARGS__);             \
        exit(1);                                         \
    }

#define RELATIVE_DIFFERENCE_FACTOR 0.2
// Inspiration: https://stackoverflow.com/a/67843030/11023871
int check_relative_difference(double number_1, double number_2) {
    if (number_1 == number_2 && number_1 == 0) {
        return 1;
    }
    double greater_magnitude =
        std::fmax(std::abs(number_1), std::abs(number_2));

    return std::abs(number_1 - number_2) <
           RELATIVE_DIFFERENCE_FACTOR * greater_magnitude;
}

void fill_with_random(double* A, int size) {
    for (int i = 0; i < size; i++) {
        A[i] = rand() % 1000;
    }
}

void fill_with_same_values(double* A, arma::mat* B, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            (*B)(i, j) = A[i * cols + j];
        }
    }
}

void fill_with_zero(double* A, int size) {
    for (int i = 0; i < size; i++) {
        A[i] = 0;
    }
}

void fill_with_random(arma::mat* A, int size) {
    for (int i = 0; i < size; i++) {
        A[i] = rand() % 1000;
    }
}

void print_matrix(double* A, int m, int n) {
    printf("----------\n");
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%lf ", A[i * n + j]);
        }
        putchar('\n');
    }
    printf("----------\n");
}

void print_matrix(arma::mat& A, int m, int n) {
    printf("----------\n");
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%lf ", A(i, j));
        }
        putchar('\n');
    }
    printf("----------\n");
}

vector<string> get_file_names() {
    // Fast tests
    vector<string> file_names;

    // Add cube tests
    for (int i = 1; i <= 15; i++) {
        file_names.push_back("./examples/cubes/cube_" + to_string(i));
    }

    // Add simplex tests
    for (int i = 1; i <= 15; i++) {
        file_names.push_back("./examples/simplex/simplex_" + to_string(i));
    }

    vector<string> file_names_robust_test = {
            "./examples/cc_8_10",
            "./examples/cross_7",
            "./examples/cross_9",
            "./examples/cube_10_2",
            "./examples/cube_10",
            "./examples/cube_14_2",
            "./examples/cube_14",
            "./examples/cube_15",
            "./examples/cube_20",
            "./examples/ex_1",
            "./examples/rect_3",
            "./examples/rh_1",
            "./examples/rh_2",
            "./examples/rh_20_40",
    };

    copy(file_names_robust_test.begin(),file_names_robust_test.end(), back_inserter(file_names));

    // If we want bigger tests
    /*vector<string> file_names = {
            "./examples/cc_8_10",
            "./examples/cc_8_11",
            "./examples/cross_7",
            "./examples/cross_9",
            "./examples/cross_13", // takes too long
            "./examples/cube_10_2",
            "./examples/cube_14_2",
            "./examples/ex_1",
            "./examples/ex_2",
            "./examples/fm_6",
            "./examples/rect_3",
            "./examples/rh_1",
            "./examples/rh_2",
            "./examples/rh_3",
            "./examples/rh_4",
            "./examples/rh_20_40",
            "./examples/rh_30_60",
            "./examples/rh_40_80",
    };*/

    return file_names;
}

double rand_double_test(double upper_bound) { return 1.0; }
int rand_int_test(int upper_bound) { return 1; }

#endif  // TEST_FRAMEWORK_H_