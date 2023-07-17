#ifndef POLYTOPE_H_
#define POLYTOPE_H_

#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "glpk.h"
#include "math.h"
#include "../openblas/linalg.h"

#define DOUBLE_EPS ((double)1e-6)
#define PI 3.1415926536

typedef struct polytope {
    size_t m, n;

    double* A;
    double* b;

    // TODO: look what this init_radius is
    double init_radius;
    // Upper bound on dimension to compute volume in
    int dimension_limit;
    // Ball radiuses in dimension_limit size
    double* ball_radiuses;
    // Sampling starting points
    double* sampling_points_start;
    // They are using substitute B and Ai for when the values get updated in
    // preprocessing. They are pointers to vecs and matrices for each possible
    // dimension
    double** A_dimensions;
    double** b_dimensions;
    // Reciprocal of beta, beta-cut
    double beta_r;

    // approximating volume variables
    double volume;
    double determinant;

    bool debug_msg;

} polytope;

void init_polytope(polytope** p, size_t m, size_t n);
void cleanup(polytope* p);
void check_halfplanes(polytope* p);

// the preprocess() from the paper
void preprocess(polytope* p);

// coef - number of balls
double estimate_volume(polytope* p, int coef);

// coef - number of balls
// Wrapper
double estimate_volume_with_rand_functions(polytope* p, int coef, double (*rand_double)(double),
                       int (*rand_int)(int));

// walk called by wrapper and tests
double walk_with_rand_functions(polytope* p, int k,
                                double (*rand_double)(double),
                                int (*rand_int)(int));
// get a new random point to sample. (wrapper)
double walk(polytope* p, int k);

// TA recommends using https://prng.di.unimi.it/ <- very fast, easy to use with
// vector instructions
double rand_double(double upper_bound);
int rand_int(int upper_bound);

void read_polytope(polytope** p, FILE* fin);
void print_polytope(polytope* p);

#endif  // POLYTOPE_H_