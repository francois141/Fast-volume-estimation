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
#include "linalg_optimized.h"
#include "math.h"

#define DOUBLE_EPS ((double)1e-6)
#define PI 3.1415926536

typedef struct polytope {
    size_t m, n;

    double *A;
    double *b;

    // TODO: look what this init_radius is
    double init_radius;
    // Upper bound on dimension to compute volume in
    int dimension_limit;
    // Ball radiuses in dimension_limit size
    double *ball_radiuses;
    // Sampling starting points
    double *sampling_points_start;
    // They are using substitute B and Ai for when the values get updated in
    // preprocessing. They are pointers to vecs and matrices for each possible
    // dimension
    double **A_dimensions;
    double **b_dimensions;
    // Reciprocal of beta, beta-cut
    double beta_r;

    // approximating volume variables
    double volume;
    double determinant;

    bool debug_msg;

    // Optimization - precompute c's
    double c1, c2, c3, c4;
} polytope;

void init_polytope(polytope **pp, size_t m, size_t n);

void cleanup(polytope *p);

void check_halfplanes(polytope *p);

bool generate_init_ellipsoid(polytope *p, double *r2, double *origin_arr);
void preprocess(polytope *p);

double uball_volume(int n);
double estimate_volume(polytope *p, int coef);
double estimate_volume_with_rand_functions(polytope *p, int coef,
                                           double (*rand_double)(double),
                                           int (*rand_int)(int));
double walk_with_rand_functions(polytope *p, int k,
                                double (*rand_double)(double),
                                int (*rand_int)(int));
double walk(polytope *p, int k);
double rand_double(double upper_bound);
int rand_int(int upper_bound);
void read_polytope(polytope **p, FILE *fin);
void print_polytope(polytope *p);

/*
 * assuming p is pointing to NULL
 * doing this because its bad practice to instantiate empty object polygon p. it
 * does a default initialization of the pointers inside the structure, which
 * causes memory leaks it is better to declare polygon *p = NULL, and then let
 * init_polytope allocate it
 */
void init_polytope(polytope **pp, size_t m, size_t n) {
    assert(*pp == NULL);

    polytope *p = (polytope *)malloc(sizeof(polytope));

    srand((int)time(
        0));  // Probably we'll move this (but both main and test need it)

    p->m = m;
    p->n = n;
    p->A = (double *)calloc(m * n, sizeof(double));
    p->b = (double *)calloc(m, sizeof(double));
    p->debug_msg = true;

    // Init radius
    p->init_radius = 2 * p->n;

    // Apply logarithm base conversion
    // They add a (+2) for dimensions of vectors, but we could have a nicer way
    p->dimension_limit = (int)ceil(p->n * log(p->init_radius) / log(2.0)) + 1;
    p->ball_radiuses = (double *)calloc(p->dimension_limit, sizeof(double));
    for (int i = 0; i < p->dimension_limit; ++i) {
        // Formula given on page 3.
        p->ball_radiuses[i] = pow(2.0, (2.0 * i) / (double)p->n);
    }

    // This is changed during walks and preprocessing
    p->sampling_points_start = (double *)calloc(p->n, sizeof(double));

    // Vector and matrix for each dimension
    p->A_dimensions = (double **)calloc(n, sizeof(double *));
    p->b_dimensions = (double **)calloc(n, sizeof(double *));
    // Initialize beta variable
    p->beta_r = 2 * p->n;
    p->volume = 0;
    p->determinant = 0;
    for (int i = 0; i < n; ++i) {
        p->A_dimensions[i] = (double *)calloc(m * n, sizeof(double));
        p->b_dimensions[i] = (double *)calloc(m, sizeof(double));
    }
    p->c1 = (2 * pow(p->n, 2) + pow(1 - p->n / p->beta_r, 2)) *
            (1 - 1.0 / pow(p->beta_r, 2)) / (2 * pow(p->n, 2) - 2);
    p->c2 = (1 - p->n / p->beta_r) / (p->n + 1);
    p->c3 = p->beta_r * p->beta_r;
    p->c4 = 2 * p->c2 / (1 - 1.0 / p->beta_r);
    *pp = p;
    return;
}

void cleanup(polytope *p) {
    free(p->A);
    free(p->b);

    free(p->ball_radiuses);
    free(p->sampling_points_start);

    for (int i = 0; i < p->n; ++i) {
        free(p->A_dimensions[i]);
        free(p->b_dimensions[i]);
    }

    free(p->A_dimensions);
    free(p->b_dimensions);

    free(p);
    return;
}

void check_halfplanes(polytope *p) {
    // init GLPK
    glp_prob *lp;

    // Optimization - get m and n
    size_t m = p->m;
    size_t n = p->n;

    lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MAX);
    glp_add_rows(lp, m);
    glp_add_cols(lp, n);

    // disable msg output
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_ERR;

    // load constraints
    int *ind = (int *)malloc(sizeof(int) * (n + 1));
    double *val = (double *)malloc(sizeof(double) * (n + 1));

    for (int i = 1; i < m + 1; i++) {
        for (int j = 1; j < n + 1; j++) {
            ind[j] = j;
            val[j] = p->A[(i - 1) * n + j - 1];
        }
        glp_set_mat_row(lp, i, n, ind, val);
        glp_set_row_bnds(lp, i, GLP_UP, 0, p->b[i - 1]);
    }

    free(ind);
    free(val);

    for (int i = 1; i < n + 1; i++) {
        glp_set_col_bnds(lp, i, GLP_FR, 0, 0);
    }

    // feasiblity check
    int num[2];
    for (int i = 1; i < m + 1;) {
        glp_set_row_bnds(lp, i, GLP_LO, p->b[i - 1] + 0.00001, 0);
        glp_set_obj_coef(lp, 1, 1);
        for (int j = 1; j < n; j++) {
            glp_set_obj_coef(lp, j + 1, 0);
        }
        glp_simplex(lp, &parm);
        if (glp_get_status(lp) == GLP_NOFEAS) {
            num[1] = i;
            glp_del_rows(lp, 1, num);
            // A.shed_row(i - 1);
            remove_row(m, n, p->A, i - 1);
            // b.shed_row(i - 1);
            remove_cell(m, p->b, i - 1);
            m--;
        } else {
            glp_set_row_bnds(lp, i, GLP_UP, 0, p->b[i - 1]);
            i++;
        }
    }

    // remove wasted memory left cause by remove_row and remove_cell
    p->A = (double *)realloc(p->A, m * n * sizeof(double));
    p->b = (double *)realloc(p->b, m * sizeof(double));

    // cout << "Hyperplanes Left: " << m << endl;
    if (p->debug_msg) printf("Hyperplanes Left: %lu\n", m);
    p->m = m;
    p->n = n;
    glp_delete_prob(lp);
    return;
}

bool generate_init_ellipsoid(polytope *p, double *r2, double *origin_arr) {
    size_t m = p->m;
    size_t n = p->n;
    *r2 = 0;
    memset(origin_arr, 0, n * sizeof(double));

    // init GLPK
    glp_prob *lp;
    lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MAX);
    glp_add_rows(lp, m);
    glp_add_cols(lp, n);

    // disable msg output
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_ERR;

    // load constraints
    int *ind = (int *)malloc((n + 1) * sizeof(int));
    double *val = (double *)malloc((n + 1) * sizeof(double));

    for (int i = 1; i < m + 1; i++) {
        for (int j = 1; j < n + 1; j++) {
            ind[j] = j;
            val[j] = p->A[(i - 1) * n + j - 1];
        }
        glp_set_mat_row(lp, i, n, ind, val);
        glp_set_row_bnds(lp, i, GLP_UP, 0, p->b[i - 1]);
    }

    free(ind);
    free(val);

    for (int i = 1; i < n + 1; i++) {
        glp_set_col_bnds(lp, i, GLP_FR, 0, 0);
    }

    // get bounds
    for (int i = 0; i < n; i++) {
        double max, min;
        for (int j = 0; j < n; j++) {
            glp_set_obj_coef(lp, j + 1, 0);
        }

        glp_set_obj_coef(lp, i + 1, 1);
        glp_simplex(lp, &parm);

        max = glp_get_obj_val(lp);

        for (int j = 0; j < n; j++) {
            origin_arr[j] += glp_get_col_prim(lp, j + 1);
        }

        glp_set_obj_coef(lp, i + 1, -1);
        glp_simplex(lp, &parm);
        min = -glp_get_obj_val(lp);

        for (int j = 0; j < n; j++) {
            origin_arr[j] += glp_get_col_prim(lp, j + 1);
        }

        *r2 += (max - min) * (max - min);
    }

    for (int i = 0; i < n; i++) {
        origin_arr[i] /= 2 * n;
    }

    glp_delete_prob(lp);

    return true;
}

// the preprocess() from the paper
void preprocess(polytope *p) {
    // Delete redundant hyperplanes
    check_halfplanes(p);

    // Optimization - get c's from p
    double c1, c2, c3, c4;
    c1 = p->c1;
    c2 = p->c2;
    c3 = p->c3;
    c4 = p->c4;
    // Optimization - get n
    size_t n, m;
    n = p->n;
    m = p->m;

    // Init ellipsoid (R2I, 0), Transformation = R2I, origin = 0.
    double *origin_arr = (double *)malloc(n * sizeof(double));
    double *r2 = (double *)malloc(sizeof(double));

    generate_init_ellipsoid(p, r2, origin_arr);

    double *T = (double *)calloc(n * n, sizeof(double));

    for (int i = 0; i < n; i++) {
        T[i * n + i] = *r2;
    }

    double *distance = (double *)calloc(m, sizeof(double));

    double *tm = (double *)calloc(m, sizeof(double));

    int counter = 0;

    double *temp1 = NULL;
    double *temp2 = NULL;
    while (++counter > 0) {
        int i;

        // distance = b - A * ori;

        double *temp1 = (double *)malloc(sizeof(double) * m);

        matrix_multiplication_wrapper(m, 1, n, 1, 0, p->A, origin_arr, temp1);
        matrix_matrix_sub(m, 1, p->b, temp1, distance);

        free(temp1);

        for (i = 0; i < m; i++) {
            if (distance[i] < 0) {
                // tm(i) = as_scalar(A.row(i) * T * A.row(i).t());
                // the transpose is stored the same way in memory if one
                // dimension is 1

                temp1 = (double *)malloc(sizeof(double) * n * 1);
                matrix_multiplication_wrapper(1, n, n, 1, 0, (p->A + n * i), T,
                                              temp1);

                temp2 = (double *)malloc(sizeof(double) * 1 * 1);
                matrix_multiplication_wrapper(1, 1, n, 1, 0, temp1,
                                              (p->A + n * i), temp2);

                tm[i] = temp2[0];

                free(temp1);
                free(temp2);

                break;
            }
        }

        if (i == m) {
            // check if small ellipsoid contained in polytope
            for (i = 0; i < m; i++) {
                // tm(i) = as_scalar(A.row(i) * T * A.row(i).t());
                // the transpose is stored the same way in memory if one
                // dimension is 1

                temp1 = (double *)malloc(sizeof(double) * n * 1);
                matrix_multiplication_wrapper(1, n, n, 1, 0, (p->A + n * i), T,
                                              temp1);

                temp2 = (double *)malloc(sizeof(double) * 1 * 1);
                matrix_multiplication_wrapper(1, 1, n, 1, 0, temp1,
                                              (p->A + n * i), temp2);

                tm[i] = temp2[0];

                free(temp1);
                free(temp2);

                if (c3 * distance[i] * distance[i] - tm[i] < 0) break;
            }
        }

        // terminate if E satisfies two criterions
        if (i == m) break;

        temp1 = (double *)malloc(sizeof(double) * n * 1);
        temp2 = (double *)malloc(sizeof(double) * n * 1);

        // vec t = T * A.row(i).t() / sqrt(tm(i));      temp1 = t
        matrix_multiplication_wrapper(n, 1, n, 1, 0, T, (p->A + n * i), temp1);
        // matrix_multiply_constant_inplace(n, 1, temp1, 1.0 / sqrt(tm[i]));
        matrix_divide_constant_inplace(n, 1, temp1, sqrt(tm[i]));

        // ori = ori - t * c2;      temp1 = t
        matrix_multiply_constant(n, 1, temp1, c2, temp2);
        matrix_matrix_sub_inplace(n, 1, origin_arr, temp2);

        free(temp2);
        temp2 = (double *)malloc(sizeof(double) * n * n * 1);

        // T = c1 * (T - c4 * t * t.t()); temp1 = t
        matrix_multiplication_wrapper(n, n, 1, 1, 0, temp1, temp1, temp2);
        matrix_multiply_constant_inplace(n, n, temp2, c4);
        matrix_matrix_sub_inplace(n, n, T, temp2);
        matrix_multiply_constant_inplace(n, n, T, c1);

        free(temp1);
        free(temp2);
    }

    double *Trans = (double *)malloc(n * n * sizeof(double));
    memcpy(Trans, T, n * n * sizeof(double));

    int return_code = cholesky_decomposition(n, Trans);
    assert(return_code == 0);

    // cout << Trans << endl;

    // b = beta_r * (b - A * ori);
    temp1 = (double *)malloc(m * 1 * sizeof(double));
    matrix_multiplication_wrapper(m, 1, n, 1, 0, p->A, origin_arr, temp1);
    matrix_matrix_sub_inplace(m, 1, p->b, temp1);
    matrix_multiply_constant_inplace(m, 1, p->b, p->beta_r);
    free(temp1);

    // A = A * Trans.t();
    // TODO, change the flags passed to cblas instead of actually computing the
    // transpose
    temp1 = (double *)malloc(n * n * sizeof(double));
    temp2 = (double *)malloc(m * n * sizeof(double));

    matrix_transpose(n, n, Trans, temp1);
    matrix_multiplication_wrapper(m, n, n, 1, 0, p->A, temp1, temp2);

    double *aux = temp2;
    temp2 = p->A;
    p->A = aux;
    // swap(p->A, temp2); // save time by not doing memcpy, but maybe check if
    // its the expected behaviour

    free(temp1);
    free(temp2);

    double *all_ones = (double *)malloc(1 * n * sizeof(double));  // 1 x N
    for (int i = 0; i < n; ++i) {
        all_ones[i] = 1;
    }

    for (int i = 0; i < n; i++) {
        double *A_divider =
            (double *)malloc(m * 1 * sizeof(double));  // M x 1 (col vector)

        for (int j = 0; j < m; ++j) {
            A_divider[j] = p->A[j * n + i];
        }
        // Divide each element of b with divider
        for (int j = 0; j < m; ++j) {
            // if (isfinite(p->b[j]) && isfinite(A_divider[j])){
            //     if (fabs(p->b[j]) < DOUBLE_EPS && fabs(A_divider[j]) <
            //     DOUBLE_EPS)
            //         p->b_dimensions[i][j] = NAN;
            //     else if (fabs(A_divider[j]) < DOUBLE_EPS)
            //         p->b_dimensions[i][j] = INFINITY;
            //     else
            p->b_dimensions[i][j] = p->b[j] / A_divider[j];
            // }
            // else
            //     p->b_dimensions[i][j] = NAN;
        }

        double *A_divider_matrix =
            (double *)malloc(m * n * sizeof(double));  // M x N
        matrix_multiplication_wrapper(m, n, 1, 1, 0, A_divider, all_ones,
                                      A_divider_matrix);

        for (int j = 0; j < m; ++j) {
            for (int k = 0; k < n; ++k) {
                // if (isfinite(p->A[j * n + k]) &&
                // isfinite(A_divider_matrix[j * n + k])){
                //     if (fabs(p->A[j * n + k]) < DOUBLE_EPS &&
                //     fabs(A_divider_matrix[j * n + k]) < DOUBLE_EPS)
                //         p->A_dimensions[i][j * n + k] = NAN;
                //     else if (fabs(A_divider_matrix[j * n + k]) <
                //     DOUBLE_EPS)
                //         p->A_dimensions[i][j * n + k] = INFINITY;
                //     else
                p->A_dimensions[i][j * n + k] =
                    p->A[j * n + k] / A_divider_matrix[j * n + k];
                // }
                // else
                //     p->A_dimensions[i][j * n + k] = NAN;
            }
        }
        free(A_divider);
        free(A_divider_matrix);
    }
    free(all_ones);

    // determinant = det(Trans) / pow(beta_r, n);
    // Optimization - p determinant
    double det = 1.0;
    for (int i = 0; i < n; i++) {
        det *= Trans[i * n + i];
    }

    det /= pow(p->beta_r, n);
    p->determinant = det;
    free(T);
    free(Trans);
    free(r2);
    free(origin_arr);
    free(distance);
    free(tm);

    return;
}

double uball_volume(int n) {
    double vol = 1;
    if (n % 2 == 1) {
        int k = (n - 1) / 2;
        vol *= pow(2, n);
        for (int i = 1; i < k + 1; i++) {
            vol *= PI * i;
        }
        for (int i = 1; i < n + 1; i++) {
            vol /= i;
        }
    } else {
        int k = n / 2;
        for (int i = 1; i < k + 1; i++) {
            vol *= PI / i;
        }
    }
    return vol;
}

// coef - number of balls
double estimate_volume(polytope *p, int coef) {
    return estimate_volume_with_rand_functions(p, coef, &rand_double,
                                               &rand_int);
}
double estimate_volume_with_rand_functions(polytope *p, int coef,
                                           double (*rand_double)(double),
                                           int (*rand_int)(int)) {
    int k, i, maxk, maxl = 0;
    const long step_size = coef * p->dimension_limit;
    double *alpha = (double *)calloc(p->dimension_limit, sizeof(double));
    long *volK = (long *)calloc(p->dimension_limit, sizeof(long));

    memset(p->sampling_points_start, 0, p->n * sizeof(double));

    for (k = p->dimension_limit - 2; k >= 0; k--) {
        for (i = volK[k + 1]; i < step_size; i++) {
            double m = walk_with_rand_functions(p, k, rand_double, rand_int);
            if (m < p->ball_radiuses[0]) {
                volK[0]++;
            } else if (m < p->ball_radiuses[k]) {
                volK[(int)trunc(p->n * log(m) / (log((double)2) * 2)) + 1]++;
            }
        }
        for (i = 0; i < k; i++) {
            volK[k] += volK[i];
        }
        if (volK[k] < step_size) {
            alpha[k] = (double)(step_size) / volK[k];
            matrix_divide_constant_inplace(1, p->n, p->sampling_points_start,
                                           pow((double)2, (double)1.0 / p->n));
        } else {
            alpha[k] = 1;
        }
    }

    p->volume = uball_volume(p->n) * p->determinant;

    for (i = 0; alpha[i] > 1 && i < p->dimension_limit - 1; i++) {
        p->volume *= alpha[i];
    }

    free(alpha);
    free(volK);

    return p->volume;
}

double walk_with_rand_functions(polytope *p, int k,
                                double (*rand_double)(double),
                                int (*rand_int)(int)) {
    double current_radius;
    double bound_max;
    double bound_min;
    double coefficient = 0.0;
    // Direction bounded by number of rows (dimensions)
    int random_direction = rand_int(p->n);

    for (int i = 0; i < p->n; ++i) {
        if (i != random_direction) {
            coefficient +=
                p->sampling_points_start[i] * p->sampling_points_start[i];
        }
    }
    // Get current radius
    current_radius = sqrt(p->ball_radiuses[k + 1] - coefficient);

    // Set maximum range
    bound_max = current_radius - p->sampling_points_start[random_direction];
    bound_min = -current_radius - p->sampling_points_start[random_direction];

    double *bounds = (double *)calloc(p->m, sizeof(double));
    // Matrix vector multiplication. Needed?
    for (int i = 0; i < p->m; ++i) {
        double prod = 0;
        for (int j = 0; j < p->n; ++j) {
            prod += p->A_dimensions[random_direction][i * p->n + j] *
                    p->sampling_points_start[j];
        }
        bounds[i] = p->b_dimensions[random_direction][i] - prod;
    }
    for (int i = 0; i < p->m; ++i) {
        if (p->A[i * p->n + random_direction] > 0) {
            if (bounds[i] < bound_max) {
                bound_max = bounds[i];
            }
        } else if (p->A[i * p->n + random_direction] < 0) {
            if (bounds[i] > bound_min) {
                bound_min = bounds[i];
            }
        }
    }

    double step = p->sampling_points_start[random_direction] +
                  rand_double(bound_max - bound_min) + bound_min;
    p->sampling_points_start[random_direction] = step;

    free(bounds);

    return coefficient + step * step;
}

// get a new random point to sample.
double walk(polytope *p, int k) {
    return walk_with_rand_functions(p, k, &rand_double, &rand_int);
}

// TA recommends using https://prng.di.unimi.it/ <- very fast, easy to use with
// vector instructions
double rand_double(double upper_bound) {
    double adjusted = (double)rand() / RAND_MAX;
    return __DBL_MIN__ + adjusted * (upper_bound - __DBL_MIN__);
}

int rand_int(int upper_bound) {
    // int adjusted = rand() / upper_bound;
    // double adjusted = (double)rand() / RAND_MAX;
    // return (0 + adjusted * (upper_bound - 0));
    return rand() % upper_bound;
}

void read_polytope(polytope **p, FILE *fin) {
    polytope *input = NULL;
    size_t m, n;
    int nb_returned_arguments = fscanf(fin, "%lu %lu", &m, &n);
    assert(nb_returned_arguments == 2);
    init_polytope(&input, m, n);
    for (int i = 0; i < input->m; i++) {
        for (int j = 0; j < input->n; j++) {
            nb_returned_arguments =
                fscanf(fin, "%lf", &input->A[i * input->n + j]);
            assert(nb_returned_arguments == 1);
        }
        nb_returned_arguments = fscanf(fin, "%lf", &input->b[i]);
        assert(nb_returned_arguments == 1);
    }

    *p = input;
}

void print_polytope(polytope *p) {
    printf("A:\n");
    for (int i = 0; i < p->m; i++) {
        for (int j = 0; j < p->n; j++) {
            // TODO: nice formatting with tabs
            printf("%lf ", p->A[i * p->n + j]);
        }
        printf("\n");
    }

    printf("b:\n");
    for (int i = 0; i < p->m; i++) {
        printf("%lf ", p->b[i]);
    }
    printf("\n");
}

#endif  // POLYTOPE_H_
