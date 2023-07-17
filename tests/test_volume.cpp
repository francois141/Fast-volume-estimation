#include <stdio.h>

#include "./original_version/vol.h"
#include "test_framework.h"
#include "test_polytope_init_helpers.h"
// extern "C" {
#include "linalg_optimized.h"
#include "polytope_optimized.h"
// }

// @DEPRECATED: Might be useful in the future
int test_volume_estimation(string example_name) {
    size_t m;
    size_t n;
    polytope* p = NULL;
    vol::polytope_orig p_orig = initialize_polytopes(example_name, &m, &n, p);
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

    assert_same(p->determinant, p_orig.determinant, "");
    // "C++ polytope determinant: %lf. Our polytope determinant: %lf\n",
    // p->determinant, p_orig.determinant);

    double res_orig = p_orig.EstimateVol(1600, &vol::polytope_orig::randd_test,
                                         &vol::polytope_orig::randi_test);
    double res_ours = estimate_volume_with_rand_functions(
        p, 1600, &rand_double_test, &rand_int_test);
    assert_same(res_orig, res_ours,
                "C++ polytope volume estimate: %lf. Our polytope volume "
                "estimate: %lf\n",
                res_orig, res_ours);
    cleanup(p);

    return 0;
}

int test_volume_estimation_with_random(string example_name) {
    size_t m;
    size_t n;
    polytope* p = NULL;
    vol::polytope_orig p_orig = initialize_polytopes(example_name, &m, &n, p);
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

    // for (int i = 0; i < n; ++i) {
    //     for (int j = 0; j < m; ++j) {
    //         for (int k = 0; k < n; ++k) {
    //             assert_same(p_orig.Ai[i](j, k), p->A_dimensions[i][j * n +
    //             k],
    //                         "");
    //             // "Result original A_%d: %lf. Result ours: %lf\n",
    //             // i, p_orig.Ai[i](j, k), p->A_dimensions[i][j * n + k]);
    //         }
    //         assert_same(p_orig.B[i](j), p->b_dimensions[i][j], "");
    //         // "Result original B: %lf. Result ours: %lf\n",
    //         // p_orig.B[i](j), p->b_dimensions[i][j]);
    //     }
    // }
    assert_same(p->determinant, p_orig.determinant, "");
    // "C++ polytope determinant: %lf. Our polytope determinant: %lf\n",
    // p->determinant, p_orig.determinant);

    double res_orig = p_orig.EstimateVol(1600);
    double res_ours = estimate_volume(p, 1600);
    fflush(stdout);
    assert_same(res_orig, res_ours,
                "C++ polytope volume estimate: %lf. Our polytope volume "
                "estimate: %lf\n",
                res_orig, res_ours);
    cleanup(p);

    return 0;
}

int main(void) {
    srand((unsigned)time(NULL));

    vector<string> file_names = get_file_names();

    for (int i = 0; i < (int)file_names.size(); ++i) {
        // Make color red for any failure
        check_test(test_volume_estimation_with_random(file_names[i]),
                   "Test %s succeeded\n", "Test %s failed (check messages)\n",
                   file_names[i].c_str());
    }

    return 0;
}