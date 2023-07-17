//
// Created by francois on 29.04.23.
//
#include <stdio.h>

#include "../tests/original_version/vol.h"
#include "../tests/test_framework.h"
// extern "C" {
// #include "../openblas/linalg.h"
#include "linalg_optimized.h"
// #include "../polytope.h"
#include "polytope_optimized.h"
// }

#include "tsc_x86.h"

// const string FILE_CYCLES_PATH = "./benchmark_infrastructure/cycles_cubes_ours.txt";
string FILE_CYCLES_PATH = "./benchmark_infrastructure/files_to_plot/";
// const string FILE_SIZES_PATH = "./benchmark_infrastructure/size_cubes_ours.txt";
// string FILE_SIZES_PATH = "./benchmark_infrastructure/size_simplex_ours.txt";
const unsigned int nb_iterations = 3;

vector<string> get_plot_filenames() {
    vector<string> file_names;
    for (int i = 1; i <= 50; i += 2) {
    // for (int i = 5; i <= 50; i += 5) {
        file_names.push_back("./benchmark_infrastructure/cubes/cube_" +
        // file_names.push_back("./examples/simplex/simplex_" +
                             to_string(i));
    }
    return file_names;
}

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
}

long long int measure_cycles(string example_name) {
    size_t m, n;

    double* vecb = NULL;
    double* matA = NULL;

    read_test_case(example_name, &m, &n, &vecb, &matA);

    polytope* p = NULL;
    init_polytope(&p, m, n);

    p->debug_msg = false;

    for (int i = 0; i < m; i++) {
        p->b[i] = vecb[i];

        for (int j = 0; j < n; j++) {
            p->A[i * n + j] = matA[i * n + j];
        }
    }

    myInt64 cycles;
    myInt64 start;

    start = start_tsc();

    preprocess(p);
    estimate_volume(p, 1600);

    cycles = stop_tsc(start);

    std::cout << example_name << " took " << cycles << " cycles to execute."
              << std::endl;

    free(matA);
    free(vecb);
    cleanup(p);

    return cycles;
}

void output_file(vector<long long int>& cycles_array,
                 vector<string> file_names) {
    ofstream outFileCycles;
    outFileCycles.open(FILE_CYCLES_PATH, ios::out | ios::trunc);

    // ofstream outFileSize;
    // outFileSize.open(FILE_SIZES_PATH, ios::out | ios::trunc);

    for (int i = 0; i < (int)cycles_array.size(); ++i) {
        outFileCycles << cycles_array[i] << std::endl;
        // outFileSize << file_names[i].substr(38) << std::endl;
        // outFileSize << file_names[i].substr(27) << std::endl;
    }
}

int main(int argc, char** argv) {
    srand((unsigned)time(NULL));

    vector<string> file_names = get_plot_filenames();
    assert(argc == 2);

    // string aux = ""
    FILE_CYCLES_PATH += string("cycles_") + string(argv[1]) + string(".txt");
    cout << FILE_CYCLES_PATH << endl;

    // Warm-up round
    for (int i = 0; i < (int)file_names.size();
         ++i) {  // Make color red for any failure
        measure_cycles(file_names[i]);
    }

    vector<long long int> cycles_array(file_names.size(), 0);

    for (int iteration = 0; iteration < nb_iterations; iteration++) {
        for (int i = 0; i < (int)file_names.size();
             ++i) {  // Make color red for any failure
            long long int cycles = measure_cycles(file_names[i]);
            cycles_array[i] += cycles / nb_iterations;
        }
    }

    output_file(cycles_array, file_names);

    return 0;
}