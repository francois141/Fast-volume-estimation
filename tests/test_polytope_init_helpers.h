#ifndef TEST_POLYTOPE_INIT_HELPERS_H_
#define TEST_POLYTOPE_INIT_HELPERS_H_
#include <cassert>

#include "./original_version/vol.h"

#include "linalg_optimized.h"
#include "polytope_optimized.h"

void read_test_case(string file_name, size_t* m, size_t* n, double** vecb_out,
                    double** matA_out) {
    FILE* fin = fopen(file_name.c_str(), "rt");
    assert(fin != NULL);

    int nb_returned_arguments = fscanf(fin, "%lu %lu", m, n);
    assert(nb_returned_arguments == 2);

    double* vecb = (double*)malloc((*m) * sizeof(double));
    double* matA = (double*)malloc((*m) * (*n) * sizeof(double));

    for (size_t i = 0; i < *m; i++) {
        nb_returned_arguments = fscanf(fin, "%lf", &vecb[i]);
        assert(nb_returned_arguments == 1);

        for (size_t j = 0; j < *n; j++) {
            nb_returned_arguments = fscanf(fin, "%lf", &matA[i * (*n) + j]);
            assert(nb_returned_arguments == 1);
        }
    }

    *vecb_out = vecb;
    *matA_out = matA;
    fclose(fin);
}

vol::polytope_orig initialize_polytopes(string& input_file, size_t* m,
                                        size_t* n, polytope*& p) {
    double* vecb = NULL;
    double* matA = NULL;

    read_test_case(input_file, m, n, &vecb, &matA);
    init_polytope(&p, *m, *n);
    vol::polytope_orig p_orig(*m, *n);

    p->debug_msg = false;
    p_orig.msg_off = true;

    for (int i = 0; i < *m; i++) {
        p_orig.vecb(vecb[i], i);
        p->b[i] = vecb[i];

        for (int j = 0; j < *n; j++) {
            p_orig.matA(matA[i * (*n) + j], i, j);
            p->A[i * (*n) + j] = matA[i * (*n) + j];
        }
    }

    free(matA);
    free(vecb);
    return p_orig;
}

#endif  // TEST_POLYTOPE_INIT_HELPERS_H_