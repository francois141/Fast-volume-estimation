#include <stdio.h>

#include "./original_version/vol.h"
#include "test_framework.h"
#include "test_polytope_init_helpers.h"
// extern "C" {
#include "linalg_optimized.h"
#include "polytope_optimized.h"
// }
#include "test_framework.h"

int test_example(string example_name) {
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

    cleanup(p);

    return 0;
}

int main(void) {
    srand((unsigned)time(NULL));

    vector<string> file_names = get_file_names();

    for (int i = 0; i < (int)file_names.size(); ++i) {
        // Make color red for any failure
        check_test(test_example(file_names[i]), "Test %s succeeded\n",
                   "Test %s failed (check messages)\n", file_names[i].c_str());
    }

    return 0;
}