#include <assert.h>
#include <stdio.h>
// extern "C" {
#include "linalg_optimized.h"
#include "polytope_optimized.h"
// }
#include "./original_version/vol.h"
#include "./original_version/walk.cpp"
#include "test_framework.h"
#include "test_polytope_init_helpers.h"

/* MODIFY THIS IF MORE TESTS ARE ADDED */
#define TESTS_ARRAY_LEN 2
#define RELATIVE_DIFFERENCE_FACTOR 0.0001

int test_init() {
    // M = 160, N = 80
    size_t m = 160;
    size_t n = 80;
    polytope* p = NULL;
    vol::polytope_orig p_orig(m, n);
    init_polytope(&p, m, n);
    assert(p != NULL);
    assert(p_orig.m == p->m && p_orig.n == p->n &&
           "Dimensions not equivalent\n");
    assert(p_orig.l == p->dimension_limit &&
           "Dimension limit not equivalent\n");
    for (int i = 0; i < p_orig.l; ++i) {
        if (!check_relative_difference(p_orig.r2[i], p->ball_radiuses[i])) {
            printf("Ball radius not close enough! Expected: %lf. Got: %lf\n",
                   p_orig.r2[i], p->ball_radiuses[i]);
            cleanup(p);
            return 1;
        }
    }
    cleanup(p);
    return 0;
}

int test_walk_ones() {
    size_t m = 160;
    size_t n = 80;
    int k = 10;
    polytope* p = NULL;
    vol::polytope_orig p_orig(m, n);
    init_polytope(&p, m, n);


    double res_new =
        walk_with_rand_functions(p, k, &rand_double_test, &rand_int_test);
    double res_orig = p_orig.walk_with_rand_functions(
        k, &vol::polytope_orig::randd_test, &vol::polytope_orig::randi_test,
        p_orig);
    if (!check_relative_difference(res_new, res_orig)) {
        printf("Results are off! Expected: %lf. Got: %lf\n", res_orig, res_new);
        cleanup(p);
        return 1;
    }
    cleanup(p);
    return 0;
}

int test_walk_with_values_set(string& file_name) {
    size_t m;
    size_t n;
    int k;
    polytope* p = NULL;
    vol::polytope_orig p_orig = initialize_polytopes(file_name, &m, &n, p);
    p_orig.Preprocess();
    preprocess(p);

    // assert_same(p->m, p_orig.m, "Result original m: %lu. Result ours: %lu\n",
    // p_orig.m, p->m);
    assert_same(p->m, p_orig.m, "");
    // assert_same(p->n, p_orig.n, "Result original n: %lu. Result ours: %lu\n",
    // p_orig.n, p->n);
    assert_same(p->n, p_orig.n, "");

    m = p->m;
    n = p->n;
    k = p->dimension_limit - 2;
    double res_new =
        walk_with_rand_functions(p, k, &rand_double_test, &rand_int_test);
    double res_orig = p_orig.walk_with_rand_functions(
        k, &vol::polytope_orig::randd_test, &vol::polytope_orig::randi_test,
        p_orig);
    if (!check_relative_difference(res_new, res_orig)) {
        printf("Results are off! Expected: %lf. Got: %lf\n", res_orig, res_new);
        cleanup(p);
        return 1;
    }
    cleanup(p);
    return 0;
}

int main(void) {
    /* Run all tests
     *  Current version simply calls the functions of the tests
     */
    int (*tests[TESTS_ARRAY_LEN])() = {&test_init, &test_walk_ones};
    for (int i = 0; i < TESTS_ARRAY_LEN; ++i) {
        check_test(tests[i](), "Test %d succeeded\n",
                   "Test %d failed (check messages)\n", i + 1);
    }

    vector<string> file_names = get_file_names();

    for (int i = 0; i < (int)file_names.size(); ++i) {
        check_test(test_walk_with_values_set(file_names[i]),
                   "Test %s succeeded\n", "Test %s failed (check messages)\n",
                   file_names[i].c_str());
    }
    return 0;
}