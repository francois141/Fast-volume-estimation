# List of optimizations

### Optimmization 1 - remove calloc calls in walk_with_rand_functions

#### Description

Call calloc only once. At each iteration the value will be freed and allocated again.

#### Code

```c
// In preprocessing
p->bounds = (double *)calloc(p->m, sizeof(double));
p->alpha = (double *)calloc(p->dimension_limit, sizeof(double));
p->volK = (long *)calloc(p->dimension_limit, sizeof(long));
```

instead of

```c
// In walk_with_rand_functions
bounds = (double *)calloc(p->m, sizeof(double));
alpha = (double *)calloc(p->dimension_limit, sizeof(double));
volK = (long *)calloc(p->dimension_limit, sizeof(long));
```

### Optimization 2 - precompute uball volume

#### Description

We can precompute the function uball_volume and just call the computed value inside of estimate_volume_with_rand_functions. The function has no side effects and always have the same input value.

#### Code

```c
p->volume = p->uball_volume * p->determinant;
```

instead of

```c
p->volume = uball_volume(p->n) * p->determinant;
```

### Optimization 3 - simplify expression

#### Description

The expression can be simplified and optimized.

#### Code

```c
volK[(int)trunc(log(m) * p->n_mult_two_log_two_inv) + 1]++;
```

instead of

```c
volK[(int)trunc(p->n * log(m) / (log((double)2) * 2)) + 1]++;
```

### Optimization 4 - explicit matrix-vector multiplication (only in preprocess)

#### Description

In `preprocess`, any matrix or vector manipulation was done using calls to our wrappers. A matrix-vector multiplication was handled using matrix-matrix multiplication, which has a higher complexity and might result in unnecessary flops. More so, a lot of the times, the operations were done through multiple calls in steps, when they could be done in one step (see example below).

#### Code

```c
for (int j = 0; j < p->n; j++) {
    for (int k = 0; k < p->n; k++) {
        T[j * p->n + k] = c1 * (T[j * p->n + k] - c4 * temp1[j] * temp1[k]); 
    }
}
```

instead of

```c
matrix_multiplication_wrapper(p->n, p->n, 1, 1, 0, temp1, temp1, temp2);
matrix_multiply_constant_inplace(p->n, p->n, temp2, c4);
matrix_matrix_sub_inplace(p->n, p->n, T, temp2);
matrix_multiply_constant_inplace(p->n, p->n, T, c1);
```
### Optimization 5 - Scalar replacement

In various places we access values through references, especially `p->n` and `p->m`. Moreover, we also update values of such references in place. We instead replace such reference usages by making use of local variables and updating references at the end.
### Optimization 6 - using the xoshiro random number generator

#### Description

Xoshiro is the random number generator recommended by the TA, from [here](https://prng.di.unimi.it/).

#### Code

```c
double rand_double(double upper_bound) {
    double adjusted = (next() >> 11) * 0x1.0p-53;
    return __DBL_MIN__ + adjusted * (upper_bound - __DBL_MIN__);
}

int rand_int(int upper_bound) {
    double adjusted = (next() >> 11) * 0x1.0p-53;
    return (0 + adjusted * (upper_bound - 0));
}
```

instead of

```c
double rand_double(double upper_bound) {
    double adjusted = (double)rand() / RAND_MAX;
    return __DBL_MIN__ + adjusted * (upper_bound - __DBL_MIN__);
}

int rand_int(int upper_bound) {
    double adjusted = (double)rand() / RAND_MAX;
    return (0 + adjusted * (upper_bound - 0));
}
```

### Optimization 7 - precompute A \* sampling_points_start product and reuse it instead of MVM in walk

#### Description

`A_dimensions[random_direction] * sampling_points_start = A * sampling_points_start / A.col(random_direction)`

#### Code

```c
for (int i = 0; i < m; ++i) {
    double prod = p->A_prod_sampling_points[i];
    prod = (p->b[i] - prod) / p->A[i * n + random_direction];
}
.
.
.
double incr = rand_double(bound_max - bound_min) + bound_min;
double step = p->sampling_points_start[random_direction] + incr;

for (int i = 0; i < m; i++){
    p->A_prod_sampling_points[i] += incr * p->A[i * n + random_direction];
}
```

instead of

```c
for (int i = 0; i < p->m; ++i) {
    double prod = 0;
    for (int j = 0; j < p->n; ++j) {
        prod += p->A_dimensions[random_direction][i * p->n + j] *
                p->sampling_points_start[j];
    }
    bounds[i] = p->b_dimensions[random_direction][i] - prod;
}
```


### Optimization 8 - keep whole p->bounds precomputed, and remove computation of coefficient in walk

#### Description

Same idea as optimization 7, but we keep the value `p->b[i] - p->A_prod_sampling_points[i]` precomputed and updating it when needed. Additionally, removed, the for loop for the coefficient.

#### Code

```c
double coefficient = p->sampling_poinst_sq_sum - p->sampling_points_start[random_direction] * p->sampling_points_start[random_direction];
...
for (int i = 0; i < m; ++i) {
    double prod = p->bounds[i] / p->A[i * n + random_direction];
    ...
}
...
for (int i = 0; i < m; i++){
    p->bounds[i] -= incr * p->A[i * n + random_direction];
}
p->sampling_poinst_sq_sum += step * step - p->sampling_points_start[random_direction] * p->sampling_points_start[random_direction];
p->sampling_points_start[random_direction] = step;
```

instead of

```c
for (int i = 0; i < p->n; ++i) {
    if (i != random_direction) {
        coefficient +=
                p->sampling_points_start[i] * p->sampling_points_start[i];
    }
}
...
for (int i = 0; i < m; ++i) {
    double prod = p->A_prod_sampling_points[i];
    prod = (p->b[i] - prod) / p->A[i * n + random_direction];
}
.
.
.
double incr = rand_double(bound_max - bound_min) + bound_min;
double step = p->sampling_points_start[random_direction] + incr;

for (int i = 0; i < m; i++){
    p->A_prod_sampling_points[i] += incr * p->A[i * n + random_direction];
}
```

### Optimization 9 - Transpose matrix and apply SIMD in the critical part of the code

We can see that this code has two problems. First, the matrix A has no spatial locality. We can optimize this by transposing the matrix. 
The second part is that we can use SIMDs to go faster. When m is not equal to 4 we can pad the vector with values that will not change the final value. 
So we can make calculations on multiples of 4 without loss of generality.

```c
__m256d bound_max_vector = _mm256_set1_pd(bound_max);
__m256d bound_min_vector = _mm256_set1_pd(bound_min);
__m256d zero_vector = _mm256_setzero_pd();

for (int i = 0; i < m_rounded; i += 4) {
    __m256d bounds_vector = _mm256_loadu_pd(p->bounds + i);
    __m256d A_transpose_vector = _mm256_loadu_pd(p->A_transpose_padded + random_direction * m_rounded + i);
    
    __m256d prod_vector = _mm256_div_pd(bounds_vector, A_transpose_vector);

    __m256d mask1 = _mm256_cmp_pd(A_transpose_vector, zero_vector, _CMP_GT_OQ);
    __m256d mask2 = _mm256_cmp_pd(A_transpose_vector, zero_vector, _CMP_LT_OQ);
    
    __m256d min_vector = _mm256_min_pd(bound_max_vector, prod_vector);
    __m256d max_vector = _mm256_max_pd(bound_min_vector, prod_vector);
    
    bound_max_vector = _mm256_blendv_pd(bound_max_vector, min_vector, mask1);
    bound_min_vector = _mm256_blendv_pd(bound_min_vector, max_vector, mask2);
}

bound_min = SIMD_REDUCE_MAX(bound_min_vector);
bound_max = SIMD_REDUCE_MIN(bound_max_vector);

double incr = rand_double(bound_max - bound_min) + bound_min;
double step = p->sampling_points_start[random_direction] + incr;

__m256d step_vector = _mm256_set1_pd(incr);
for (int i = 0; i < m_rounded; i+=4){
    // p->bounds[i] -= incr * p->A_transpose_padded[random_direction * m_rounded + i];
    __m256d bounds_vector = _mm256_loadu_pd(p->bounds + i);
    __m256d A_transpose_vector = _mm256_loadu_pd(p->A_transpose_padded + random_direction * m_rounded + i);
    _mm256_storeu_pd(p->bounds + i, _mm256_fnmadd_pd(step_vector,A_transpose_vector, bounds_vector));
}
```

instead of

```c
for (int i = 0; i < m; ++i) {
    double prod = p->bounds[i] / p->A[i * n + random_direction];
    
    if (p->A[i * p->n + random_direction] > 0) {
        if (prod < bound_max) {
            bound_max = prod;
        }
    } else if (p->A[i * p->n + random_direction] < 0) {
            if (prod > bound_min) {
                bound_min = prod;
            }
        }
    }
    
    double incr = rand_double(bound_max - bound_min) + bound_min;
    double step = p->sampling_points_start[random_direction] + incr;
    
    for (int i = 0; i < m; i++){
    p->bounds[i] -= incr * p->A[i * n + random_direction];
}
```

### Optimization 10 - loop unrolling

For the `walk_with_rand_functions` function where we make use of SIMD instructions, we can try unrolling with a factor of 4, 2 and 1. Due to sizes being multiples of 4, we can't always unroll with a factor of 4 or even 2, so we must add subsequent loops to handle these cases.

Another small optimization was to initialize the `zero_vector` once for each `estimate_volume` call.

#### Code

```c
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
```

instead of:

```c
__m256d zero_vector = _mm256_setzero_pd();

    for (int i = 0; i < m_rounded; i += 4) {
        __m256d bounds_vector = _mm256_loadu_pd(p->bounds + i);
        __m256d A_transpose_vector = _mm256_loadu_pd(
            p->A_transpose_padded + random_direction * m_rounded + i);

        // double prod = p->bounds[i] / p->A_transpose_padded[random_direction *
        // m_rounded + i];
        __m256d prod_vector = _mm256_div_pd(bounds_vector, A_transpose_vector);

        /*if (p->A_transpose_padded[random_direction * m_rounded + i] > 0) {
            bound_max = MAX(bound_max, prod);
        } else if (p->A_transpose_padded[random_direction * m_rounded + i] < 0)
        { bound_min = MAX(bound_min, prod);
        }*/
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
    for (int i = 0; i < m_rounded; i += 4) {
        // p->bounds[i] -= incr * p->A_transpose_padded[random_direction *
        // m_rounded + i];
        __m256d bounds_vector = _mm256_loadu_pd(p->bounds + i);
        __m256d A_transpose_vector = _mm256_loadu_pd(
            p->A_transpose_padded + random_direction * m_rounded + i);
        _mm256_storeu_pd(
            p->bounds + i,
            _mm256_fnmadd_pd(step_vector, A_transpose_vector, bounds_vector));
    }
```

### Optimization 11 - Apply SIMD in method estimate_volume

In the method `estimate_volume`, there are 3 loops which could be improved using loop unrolling and SIMD intrisincs. The first loop however doesn't benefit from any optimizations (loop unrolling with different accumulators OR simd), so leave unchanged.

For SIMD, we unrolled by 4. Any remaining elements are treated in a separate loop.

#### Code

```c
        if (volK[k] < step_size) {
            alpha[k] = (double)(step_size) / volK[k];

            {
                // bounds[i1] = b[i1] - (b[i1] - bounds[i1]) * ct
                int i1;
                for (i1 = 0; i1 < m-3; i1 += 4) {
                    __m256d bounds_vector = _mm256_load_pd(p->bounds + i1);
                    __m256d b_vector = _mm256_loadu_pd(p->b + i1);

                    __m256d tmp1 = _mm256_sub_pd(b_vector, bounds_vector);
                    bounds_vector = _mm256_fnmadd_pd(tmp1, pow_two_inv_n_inv_vector, b_vector);

                    _mm256_storeu_pd(p->bounds + i1, bounds_vector);
                }

                for (; i1 < m; i1++) {
                    p->bounds[i1] = p->b[i1] - (p->b[i1] - p->bounds[i1]) * p->pow_two_inv_n_inv;
                }
            }

            {
                __m256d sampling_poinst_sq_sum_vector = _mm256_set1_pd(0);
                int i1;

                for (i1 = 0; i1 < n-3; i1 += 4) {
                    __m256d sampling_points_vector = _mm256_load_pd(p->sampling_points_start + i1);
                    sampling_points_vector = _mm256_mul_pd(sampling_points_vector, pow_two_inv_n_inv_vector);

                    sampling_poinst_sq_sum_vector =
                        _mm256_fmadd_pd(sampling_points_vector, sampling_points_vector, sampling_poinst_sq_sum_vector);

                    _mm256_storeu_pd(p->sampling_points_start + i1, sampling_points_vector);
                }
                double sq_sum[4];
                _mm256_storeu_pd(&sq_sum[0], sampling_poinst_sq_sum_vector);
                p->sampling_poinst_sq_sum = sq_sum[0] + sq_sum[1] + sq_sum[2] + sq_sum[3];
                
                for (; i1 < n; i1++) {
                    p->sampling_points_start[i1] *= p->pow_two_inv_n_inv;
                    p->sampling_poinst_sq_sum += p->sampling_points_start[i1] * p->sampling_points_start[i1];
                }
            }

        } else {
            alpha[k] = 1;
        }
```

instead of

```c
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
```



### Optimization 12 - Apply SIMD in method estimate_volume with padded arrays

Similar to optim. 11.

In the method `estimate_volume`, there are 3 loops which could be improved using loop unrolling and SIMD intrisincs. The first loop however doesn't benefit from any optimizations (loop unrolling with different accumulators OR simd), so leave unchanged.

For SIMD, we unrolled by 4. The arrays are padded to next multiple of 4.

#### Code

```c
        if (volK[k] < step_size) {
            alpha[k] = (double)(step_size) / volK[k];

            {
                // bounds[i1] = b[i1] - (b[i1] - bounds[i1]) * ct
                int i1;
                for (i1 = 0; i1 < m_rounded; i1 += 4) {
                    __m256d bounds_vector = _mm256_load_pd(p->bounds + i1);
                    __m256d b_vector = _mm256_loadu_pd(p->b + i1);

                    __m256d tmp1 = _mm256_sub_pd(b_vector, bounds_vector);
                    bounds_vector = _mm256_fnmadd_pd(tmp1, pow_two_inv_n_inv_vector, b_vector);

                    _mm256_storeu_pd(p->bounds + i1, bounds_vector);
                }
            }

            {
                __m256d sampling_poinst_sq_sum_vector = _mm256_set1_pd(0);
                int i1;

                for (i1 = 0; i1 < n_rounded; i1 += 4) {
                    __m256d sampling_points_vector = _mm256_load_pd(p->sampling_points_start + i1);
                    sampling_points_vector = _mm256_mul_pd(sampling_points_vector, pow_two_inv_n_inv_vector);

                    sampling_poinst_sq_sum_vector =
                        _mm256_fmadd_pd(sampling_points_vector, sampling_points_vector, sampling_poinst_sq_sum_vector);

                    _mm256_storeu_pd(p->sampling_points_start + i1, sampling_points_vector);
                }
                double sq_sum[4];
                _mm256_storeu_pd(&sq_sum[0], sampling_poinst_sq_sum_vector);
                p->sampling_poinst_sq_sum = sq_sum[0] + sq_sum[1] + sq_sum[2] + sq_sum[3];
            }

        } else {
            alpha[k] = 1;
        }
```

instead of

```c
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
```


