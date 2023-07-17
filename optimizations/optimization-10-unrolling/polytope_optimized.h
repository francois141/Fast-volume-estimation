#ifndef POLYTOPE_H_
#define POLYTOPE_H_

#include <assert.h>
#include <immintrin.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <cfloat>

#include "../xoshiro.h"
#include "glpk.h"
#include "linalg_optimized.h"
#include "math.h"

#define ROUND_UP(N, S) ((((N) + (S)-1) / (S)) * (S))

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

#define SIMD_REDUCE_ADD(x) (x[0] + x[1] + x[2] + x[3])
#define SIMD_REDUCE_MIN(x) (MIN(MIN(x[0], x[1]), MIN(x[2], x[3])))
#define SIMD_REDUCE_MAX(x) (MAX(MAX(x[0], x[1]), MAX(x[2], x[3])))

#define DOUBLE_EPS ((double)1e-6)
#define PI 3.1415926536

typedef struct polytope {
    size_t m, n;

    double *A;
    double *b;
    double *A_transpose_padded;

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

    // Custom precomputed values
    double uball_volume;
    double two_log_two_inv;
    double pow_two_inv_n_inv;
    double *bounds;

    double sampling_poinst_sq_sum;

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
                                int (*rand_int)(int), __m256d zero_vector);
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
    s[0] = time(0);
    s[1] = time(0);
    s[2] = time(0);
    s[3] = time(0);

    p->m = m;
    p->n = n;
    p->A = (double *)calloc(m * n, sizeof(double));
    p->b = (double *)calloc(m, sizeof(double));
    p->debug_msg = true;

    // Init radius
    p->init_radius = 2 * n;

    // Apply logarithm base conversion
    // They add a (+2) for dimensions of vectors, but we could have a nicer way
    p->dimension_limit = (int)ceil(n * log(p->init_radius) / log(2.0)) + 1;
    p->ball_radiuses = (double *)calloc(p->dimension_limit, sizeof(double));
    for (int i = 0; i < p->dimension_limit; ++i) {
        // Formula given on page 3.
        p->ball_radiuses[i] = pow(2.0, (2.0 * i) / (double)n);
    }

    // This is changed during walks and preprocessing
    p->sampling_points_start = (double *)calloc(n, sizeof(double));

    // Vector and matrix for each dimension
    p->A_dimensions = (double **)calloc(n, sizeof(double *));
    p->b_dimensions = (double **)calloc(n, sizeof(double *));
    // Initialize beta variable
    p->beta_r = 2 * n;
    p->volume = 0;
    p->determinant = 0;
    for (int i = 0; i < n; ++i) {
        p->A_dimensions[i] = (double *)calloc(m * n, sizeof(double));
        p->b_dimensions[i] = (double *)calloc(m, sizeof(double));
    }
    p->c1 = (2 * pow(n, 2) + pow(1 - n / p->beta_r, 2)) *
            (1 - 1.0 / pow(p->beta_r, 2)) / (2 * pow(n, 2) - 2);
    p->c2 = (1 - n / p->beta_r) / (n + 1);
    p->c3 = p->beta_r * p->beta_r;
    p->c4 = 2 * p->c2 / (1 - 1.0 / p->beta_r);

    p->sampling_poinst_sq_sum = 0;

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
    free(p->bounds);

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
    double r2;

    generate_init_ellipsoid(p, &r2, origin_arr);

    double *T = (double *)calloc(n * n, sizeof(double));

    for (int i = 0; i < n; i++) {
        T[i * n + i] = r2;
    }

    double *distance = (double *)calloc(m, sizeof(double));

    double *tm = (double *)calloc(m, sizeof(double));

    int counter = 0;

    double *temp1 = (double *)malloc(sizeof(double) * n * 1);
    while (++counter > 0) {
        int i;

        // distance = b - A * ori;
        for (int i = 0; i < m; i++) {
            double temp = 0.0;
            for (int j = 0; j < n; j++) {
                temp += p->A[i * n + j] * origin_arr[j];
            }
            distance[i] = p->b[i] - temp;
        }

        for (i = 0; i < m; i++) {
            if (distance[i] < 0) {
                // tm(i) = as_scalar(A.row(i) * T * A.row(i).t());
                // the transpose is stored the same way in memory if one
                // dimension is 1

                double temp = 0.0;
                for (int j = 0; j < n; j++) {
                    double temp2 = 0.0;
                    for (int k = 0; k < n; k++) {
                        temp2 += p->A[i * n + k] * T[k * n + j];
                    }
                    temp += temp2 * p->A[i * n + j];
                }
                tm[i] = temp;

                break;
            }
        }

        if (i == m) {
            // check if small ellipsoid contained in polytope
            for (i = 0; i < m; i++) {
                // tm(i) = as_scalar(A.row(i) * T * A.row(i).t());
                // the transpose is stored the same way in memory if one
                // dimension is 1

                double temp = 0.0;
                for (int j = 0; j < n; j++) {
                    double temp2 = 0.0;
                    for (int k = 0; k < n; k++) {
                        temp2 += p->A[i * n + k] * T[k * n + j];
                    }
                    temp += temp2 * p->A[i * n + j];
                }
                tm[i] = temp;

                if (c3 * distance[i] * distance[i] - tm[i] < 0) break;
            }
        }

        // terminate if E satisfies two criterions
        if (i == m) break;

        // vec t = T * A.row(i).t() / sqrt(tm(i));      temp1 = t
        double aux = 1.0 / sqrt(tm[i]);
        for (int j = 0; j < n; j++) {
            double temp = 0.0;
            for (int k = 0; k < n; k++) {
                temp += T[j * n + k] * p->A[i * n + k];
            }
            // temp1[j] = temp / sqrt(tm[i]);
            temp1[j] = temp * aux;
        }

        // ori = ori - t * c2;      temp1 = t
        for (int j = 0; j < n; j++) {
            origin_arr[j] -= temp1[j] * c2;
        }

        // T = c1 * (T - c4 * t * t.t()); temp1 = t
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                T[j * n + k] = c1 * (T[j * n + k] - c4 * (temp1[j] * temp1[k]));
            }
        }
    }

    free(temp1);

    double *Trans = (double *)malloc(n * n * sizeof(double));
    memcpy(Trans, T, n * n * sizeof(double));

    int return_code = cholesky_decomposition(n, Trans);
    assert(return_code == 0);

    // b = beta_r * (b - A * ori);
    /* // Somehow it's faster with the wrappers. Leave as is.
    for (int j = 0; j < p->m; j++) {
        double temp = 0;
        for (int k = 0; k < p->n; k++) {
            temp += p->A[j * p->n + k] * origin_arr[k];
        }
        p->b[j] = p->beta_r * (p->b[j] - temp);
    }
    */
    temp1 = (double *)malloc(m * 1 * sizeof(double));
    matrix_multiplication_wrapper(m, 1, n, 1, 0, p->A, origin_arr, temp1);
    matrix_matrix_sub_inplace(m, 1, p->b, temp1);
    matrix_multiply_constant_inplace(m, 1, p->b, p->beta_r);
    free(temp1);

    // A = A * Trans.t();
    // TODO, change the flags passed to cblas instead of actually computing the
    // transpose
    temp1 = (double *)malloc(n * n * sizeof(double));
    double *temp2 = (double *)malloc(m * n * sizeof(double));

    matrix_transpose(n, n, Trans, temp1);
    matrix_multiplication_wrapper(m, n, n, 1, 0, p->A, temp1, temp2);

    double *aux = temp2;
    temp2 = p->A;
    p->A = aux;

    // Create transposed and padded matrix A
    int n_rounded = ROUND_UP(n, 4);
    int m_rounded = ROUND_UP(m, 4);
    p->A_transpose_padded =
        (double *)malloc(sizeof(double) * m_rounded * n_rounded);

    for (int i = 0; i < n_rounded; i++) {
        for (int j = 0; j < m_rounded; j++) {
            if (i >= n || j >= m) {
                p->A_transpose_padded[i * m_rounded + j] = 0;
            } else {
                p->A_transpose_padded[i * m_rounded + j] = p->A[j * n + i];
            }
        }
    }

    // swap(p->A, temp2); // save time by not doing memcpy, but maybe check if
    // its the expected behaviour

    free(temp1);
    free(temp2);

    // for (int i = 0; i < n; i++) {
    //     for (int j = 0; j < m; ++j) {
    //         double divider =  p->A[j * n + i];
    //         p->b_dimensions[i][j] = p->b[j] / divider;
    //         for (int k = 0; k < n; ++k) {
    //             p->A_dimensions[i][j * n + k] =
    //                     p->A[j * n + k] / divider;
    //         }
    //     }
    // }

    // determinant = det(Trans) / pow(beta_r, n);
    // Optimization - p determinant
    double det = 1.0;
    for (int i = 0; i < n; i++) {
        det *= Trans[i * n + i];
    }

    p->determinant = det;
    p->determinant /= pow(p->beta_r, n);

    // Custom preprocessing
    p->uball_volume = uball_volume(n);
    p->two_log_two_inv = 1.0 / (log((double)2) * 2);
    p->pow_two_inv_n_inv = 1.0 / pow((double)2, (double)1.0 / n);
    p->bounds = (double *)calloc(ROUND_UP(m, 4), sizeof(double));

    for (int i = 0; i < m; i++) {
        double acc = 0;
        for (int j = 0; j < n; j++) {
            acc += p->A[i * n + j] * p->sampling_points_start[j];
        }
        p->bounds[i] = p->b[i] - acc;
    }

    free(T);
    free(Trans);
    // free(r2);
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
    size_t m = p->m;
    size_t n = p->n;

    int k, i, maxk, maxl = 0;
    const long long step_size = 1LL * coef * p->dimension_limit;
    double *alpha = (double *)calloc(p->dimension_limit, sizeof(double));
    long *volK = (long *)calloc(p->dimension_limit, sizeof(long));
    __m256d zero_vector = _mm256_setzero_pd();
    for (k = p->dimension_limit - 2; k >= 0; k--) {
        for (i = volK[k + 1]; i < step_size; i++) {
            double m_walk = walk_with_rand_functions(p, k, rand_double,
                                                     rand_int, zero_vector);
            if (m_walk < p->ball_radiuses[0]) {
                volK[0]++;
            } else if (m_walk < p->ball_radiuses[k]) {
                volK[(int)trunc(n * log(m_walk) * p->two_log_two_inv) + 1]++;
            }
        }

        for (i = 0; i < k; i++) {
            volK[k] += volK[i];
        }
        if (volK[k] < step_size) {
            alpha[k] = (double)(step_size) / volK[k];
            for (int i1 = 0; i1 < m; i1++) {
                p->bounds[i1] = p->b[i1] - (p->b[i1] - p->bounds[i1]) *
                                               p->pow_two_inv_n_inv;
            }
            p->sampling_poinst_sq_sum = 0;
            for (int i1 = 0; i1 < n; i1++) {
                p->sampling_points_start[i1] *= p->pow_two_inv_n_inv;

                p->sampling_poinst_sq_sum +=
                    p->sampling_points_start[i1] * p->sampling_points_start[i1];
            }
        } else {
            alpha[k] = 1;
        }
    }

    p->volume = p->uball_volume * p->determinant;

    for (i = 0; alpha[i] > 1 && i < p->dimension_limit - 1; i++) {
        p->volume *= alpha[i];
    }

    free(alpha);
    free(volK);

    return p->volume;
}

double walk_with_rand_functions(polytope *p, int k,
                                double (*rand_double)(double),
                                int (*rand_int)(int), __m256d zero_vector) {
    size_t m = p->m;
    size_t n = p->n;
    size_t m_rounded = ROUND_UP(p->m, 4);
    size_t n_rounded = ROUND_UP(p->n, 4);

    double current_radius;
    double bound_max;
    double bound_min;
    // Direction bounded by number of rows (dimensions)
    int random_direction = rand_int(n);

    double coefficient = p->sampling_poinst_sq_sum -
                         p->sampling_points_start[random_direction] *
                             p->sampling_points_start[random_direction];
    // Get current radius
    current_radius = sqrt(p->ball_radiuses[k + 1] - coefficient);

    // Set maximum range
    bound_max = current_radius - p->sampling_points_start[random_direction];
    bound_min = -current_radius - p->sampling_points_start[random_direction];

    __m256d bound_max_vector = _mm256_set1_pd(bound_max);
    __m256d bound_min_vector = _mm256_set1_pd(bound_min);
    int i;
    for (i = 0; i + 16 < m_rounded; i += 16) {
        // Loads
        __m256d bounds_vector = _mm256_loadu_pd(p->bounds + i);
        __m256d A_transpose_vector = _mm256_loadu_pd(
            p->A_transpose_padded + random_direction * m_rounded + i);

        __m256d bounds_vector_2 = _mm256_loadu_pd(p->bounds + i + 4);
        __m256d A_transpose_vector_2 = _mm256_loadu_pd(
            p->A_transpose_padded + random_direction * m_rounded + i + 4);

        __m256d bounds_vector_3 = _mm256_loadu_pd(p->bounds + i + 8);
        __m256d A_transpose_vector_3 = _mm256_loadu_pd(
            p->A_transpose_padded + random_direction * m_rounded + i + 8);

        __m256d bounds_vector_4 = _mm256_loadu_pd(p->bounds + i + 12);
        __m256d A_transpose_vector_4 = _mm256_loadu_pd(
            p->A_transpose_padded + random_direction * m_rounded + i + 12);

        // Divs
        __m256d prod_vector = _mm256_div_pd(bounds_vector, A_transpose_vector);
        __m256d prod_vector_2 =
            _mm256_div_pd(bounds_vector_2, A_transpose_vector_2);
        __m256d prod_vector_3 =
            _mm256_div_pd(bounds_vector_3, A_transpose_vector_3);
        __m256d prod_vector_4 =
            _mm256_div_pd(bounds_vector_4, A_transpose_vector_4);

        // CMP
        __m256d mask_gt =
            _mm256_cmp_pd(A_transpose_vector, zero_vector, _CMP_GT_OQ);
        __m256d mask_lt =
            _mm256_cmp_pd(A_transpose_vector, zero_vector, _CMP_LT_OQ);

        __m256d mask_gt_2 =
            _mm256_cmp_pd(A_transpose_vector_2, zero_vector, _CMP_GT_OQ);
        __m256d mask_lt_2 =
            _mm256_cmp_pd(A_transpose_vector_2, zero_vector, _CMP_LT_OQ);

        __m256d mask_gt_3 =
            _mm256_cmp_pd(A_transpose_vector_3, zero_vector, _CMP_GT_OQ);
        __m256d mask_lt_3 =
            _mm256_cmp_pd(A_transpose_vector_3, zero_vector, _CMP_LT_OQ);

        __m256d mask_gt_4 =
            _mm256_cmp_pd(A_transpose_vector_4, zero_vector, _CMP_GT_OQ);
        __m256d mask_lt_4 =
            _mm256_cmp_pd(A_transpose_vector_4, zero_vector, _CMP_LT_OQ);

        // Min/Max vectors
        __m256d min_vector = _mm256_min_pd(bound_max_vector, prod_vector);
        __m256d max_vector = _mm256_max_pd(bound_min_vector, prod_vector);

        bound_max_vector =
            _mm256_blendv_pd(bound_max_vector, min_vector, mask_gt);
        bound_min_vector =
            _mm256_blendv_pd(bound_min_vector, max_vector, mask_lt);

        __m256d min_vector_2 = _mm256_min_pd(bound_max_vector, prod_vector_2);
        __m256d max_vector_2 = _mm256_max_pd(bound_min_vector, prod_vector_2);

        bound_max_vector =
            _mm256_blendv_pd(bound_max_vector, min_vector_2, mask_gt_2);
        bound_min_vector =
            _mm256_blendv_pd(bound_min_vector, max_vector_2, mask_lt_2);

        __m256d min_vector_3 = _mm256_min_pd(bound_max_vector, prod_vector_3);
        __m256d max_vector_3 = _mm256_max_pd(bound_min_vector, prod_vector_3);

        bound_max_vector =
            _mm256_blendv_pd(bound_max_vector, min_vector_3, mask_gt_3);
        bound_min_vector =
            _mm256_blendv_pd(bound_min_vector, max_vector_3, mask_lt_3);

        __m256d min_vector_4 = _mm256_min_pd(bound_max_vector, prod_vector_4);
        __m256d max_vector_4 = _mm256_max_pd(bound_min_vector, prod_vector_4);

        bound_max_vector =
            _mm256_blendv_pd(bound_max_vector, min_vector_4, mask_gt_4);
        bound_min_vector =
            _mm256_blendv_pd(bound_min_vector, max_vector_4, mask_lt_4);
    }

    for (; i + 8 < m_rounded; i += 8) {
        // Loads
        __m256d bounds_vector = _mm256_loadu_pd(p->bounds + i);
        __m256d A_transpose_vector = _mm256_loadu_pd(
            p->A_transpose_padded + random_direction * m_rounded + i);

        __m256d bounds_vector_2 = _mm256_loadu_pd(p->bounds + i + 4);
        __m256d A_transpose_vector_2 = _mm256_loadu_pd(
            p->A_transpose_padded + random_direction * m_rounded + i + 4);

        // Divs
        __m256d prod_vector = _mm256_div_pd(bounds_vector, A_transpose_vector);
        __m256d prod_vector_2 =
            _mm256_div_pd(bounds_vector_2, A_transpose_vector_2);

        // CMP
        __m256d mask_gt =
            _mm256_cmp_pd(A_transpose_vector, zero_vector, _CMP_GT_OQ);
        __m256d mask_lt =
            _mm256_cmp_pd(A_transpose_vector, zero_vector, _CMP_LT_OQ);

        __m256d mask_gt_2 =
            _mm256_cmp_pd(A_transpose_vector_2, zero_vector, _CMP_GT_OQ);
        __m256d mask_lt_2 =
            _mm256_cmp_pd(A_transpose_vector_2, zero_vector, _CMP_LT_OQ);
        // Min/Max vectors
        __m256d min_vector = _mm256_min_pd(bound_max_vector, prod_vector);
        __m256d max_vector = _mm256_max_pd(bound_min_vector, prod_vector);

        bound_max_vector =
            _mm256_blendv_pd(bound_max_vector, min_vector, mask_gt);
        bound_min_vector =
            _mm256_blendv_pd(bound_min_vector, max_vector, mask_lt);

        __m256d min_vector_2 = _mm256_min_pd(bound_max_vector, prod_vector_2);
        __m256d max_vector_2 = _mm256_max_pd(bound_min_vector, prod_vector_2);

        bound_max_vector =
            _mm256_blendv_pd(bound_max_vector, min_vector_2, mask_gt_2);
        bound_min_vector =
            _mm256_blendv_pd(bound_min_vector, max_vector_2, mask_lt_2);
    }

    for (; i < m_rounded; i += 4) {
        __m256d bounds_vector = _mm256_loadu_pd(p->bounds + i);
        __m256d A_transpose_vector = _mm256_loadu_pd(
            p->A_transpose_padded + random_direction * m_rounded + i);

        __m256d prod_vector = _mm256_div_pd(bounds_vector, A_transpose_vector);

        __m256d mask1 =
            _mm256_cmp_pd(A_transpose_vector, zero_vector, _CMP_GT_OQ);
        __m256d mask2 =
            _mm256_cmp_pd(A_transpose_vector, zero_vector, _CMP_LT_OQ);

        __m256d min_vector = _mm256_min_pd(bound_max_vector, prod_vector);
        __m256d max_vector = _mm256_max_pd(bound_min_vector, prod_vector);

        bound_max_vector =
            _mm256_blendv_pd(bound_max_vector, min_vector, mask1);
        bound_min_vector =
            _mm256_blendv_pd(bound_min_vector, max_vector, mask2);
    }

    bound_min = SIMD_REDUCE_MAX(bound_min_vector);
    bound_max = SIMD_REDUCE_MIN(bound_max_vector);

    double incr = rand_double(bound_max - bound_min) + bound_min;
    double step = p->sampling_points_start[random_direction] + incr;

    __m256d step_vector = _mm256_set1_pd(incr);
    i = 0;
    for (; i + 16 < m_rounded; i += 16) {
        // Loads
        __m256d bounds_vector = _mm256_loadu_pd(p->bounds + i);
        __m256d A_transpose_vector = _mm256_loadu_pd(
            p->A_transpose_padded + random_direction * m_rounded + i);

        __m256d bounds_vector_2 = _mm256_loadu_pd(p->bounds + i + 4);
        __m256d A_transpose_vector_2 = _mm256_loadu_pd(
            p->A_transpose_padded + random_direction * m_rounded + i + 4);

        __m256d bounds_vector_3 = _mm256_loadu_pd(p->bounds + i + 8);
        __m256d A_transpose_vector_3 = _mm256_loadu_pd(
            p->A_transpose_padded + random_direction * m_rounded + i + 8);

        __m256d bounds_vector_4 = _mm256_loadu_pd(p->bounds + i + 12);
        __m256d A_transpose_vector_4 = _mm256_loadu_pd(
            p->A_transpose_padded + random_direction * m_rounded + i + 12);

        _mm256_storeu_pd(
            p->bounds + i,
            _mm256_fnmadd_pd(step_vector, A_transpose_vector, bounds_vector));

        _mm256_storeu_pd(p->bounds + i + 4,
                         _mm256_fnmadd_pd(step_vector, A_transpose_vector_2,
                                          bounds_vector_2));

        _mm256_storeu_pd(p->bounds + i + 8,
                         _mm256_fnmadd_pd(step_vector, A_transpose_vector_3,
                                          bounds_vector_3));

        _mm256_storeu_pd(p->bounds + i + 12,
                         _mm256_fnmadd_pd(step_vector, A_transpose_vector_4,
                                          bounds_vector_4));
    }

    for (; i + 8 < m_rounded; i += 8) {
        // Loads
        __m256d bounds_vector = _mm256_loadu_pd(p->bounds + i);
        __m256d A_transpose_vector = _mm256_loadu_pd(
            p->A_transpose_padded + random_direction * m_rounded + i);

        __m256d bounds_vector_2 = _mm256_loadu_pd(p->bounds + i + 4);
        __m256d A_transpose_vector_2 = _mm256_loadu_pd(
            p->A_transpose_padded + random_direction * m_rounded + i + 4);

        _mm256_storeu_pd(
            p->bounds + i,
            _mm256_fnmadd_pd(step_vector, A_transpose_vector, bounds_vector));

        _mm256_storeu_pd(p->bounds + i + 4,
                         _mm256_fnmadd_pd(step_vector, A_transpose_vector_2,
                                          bounds_vector_2));
    }
    for (; i < m_rounded; i += 4) {
        // p->bounds[i] -= incr * p->A_transpose_padded[random_direction *
        // m_rounded + i];
        __m256d bounds_vector = _mm256_loadu_pd(p->bounds + i);
        __m256d A_transpose_vector = _mm256_loadu_pd(
            p->A_transpose_padded + random_direction * m_rounded + i);
        _mm256_storeu_pd(
            p->bounds + i,
            _mm256_fnmadd_pd(step_vector, A_transpose_vector, bounds_vector));
    }

    p->sampling_poinst_sq_sum +=
        step * step - p->sampling_points_start[random_direction] *
                          p->sampling_points_start[random_direction];
    p->sampling_points_start[random_direction] = step;

    return coefficient + step * step;
}

// get a new random point to sample.
double walk(polytope *p, int k) {
    __m256d zero_vector = _mm256_setzero_pd();
    return walk_with_rand_functions(p, k, &rand_double, &rand_int, zero_vector);
}

// TA recommends using https://prng.di.unimi.it/ <- very fast, easy to use with
// vector instructions
double rand_double(double upper_bound) {
    double adjusted = (next() >> 11) * 0x1.0p-53;
    return __DBL_MIN__ + adjusted * (upper_bound - __DBL_MIN__);
}

int rand_int(int upper_bound) {
    // double adjusted = (next() >> 11) * 0x1.0p-53;
    // return (0 + adjusted * (upper_bound - 0));
    return next() % upper_bound;
}

void read_polytope(polytope **p, FILE *fin) {
    polytope *input = NULL;
    size_t m, n;
    int nb_returned_arguments = fscanf(fin, "%lu %lu", &m, &n);
    assert(nb_returned_arguments == 2);
    init_polytope(&input, m, n);
    for (int i = 0; i < input->m; i++) {
        nb_returned_arguments = fscanf(fin, "%lf", &input->b[i]);
        assert(nb_returned_arguments == 1);
        for (int j = 0; j < input->n; j++) {
            nb_returned_arguments =
                fscanf(fin, "%lf", &input->A[i * input->n + j]);
            assert(nb_returned_arguments == 1);
        }
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
