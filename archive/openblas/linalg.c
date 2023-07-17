#include "linalg.h"

// extern void dpotrf_(const char*, int*, double*, int*, int*, size_t);

void matrix_multiplication_wrapper(int m, int n, int k, double alpha,
                                   double beta, double* A, double* B,
                                   double* C) {
    // CBLAS because of row major ordering
    return cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                       alpha, A, k, B, n, beta, C, n);
}

int cholesky_decomposition(int N, double* A) {
    const char upper = 'U';
    int error_code = 0;

    dpotrf_(&upper, &N, A, &N, &error_code, (size_t)1);

    if (error_code != 0) {
        return error_code;
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++) {
            int idx1 = N * i + j;
            int idx2 = N * j + i;
            // double aux = A[idx1];
            A[idx2] = A[idx1];
            A[idx1] = 0;
        }
    }

    return 0;
}

void matrix_multiply_constant_inplace(int m, int n, double* A, double c) {
    for (int i = 0; i < m * n; i++) A[i] *= c;
}

void matrix_divide_constant_inplace(int m, int n, double* A, double c) {
    for (int i = 0; i < m * n; i++) A[i] /= c;
}

void matrix_multiply_constant(int m, int n, double* A, double c, double* B) {
    for (int i = 0; i < m * n; i++) {
        B[i] = A[i] * c;
    }
}

void matrix_matrix_add_inplace(int m, int n, double* A, double* B) {
    for (int i = 0; i < m * n; i++) {
        A[i] += B[i];
    }
}

void matrix_matrix_add(int m, int n, double* A, double* B, double* C) {
    for (int i = 0; i < m * n; i++) {
        C[i] = A[i] + B[i];
    }
}

void matrix_matrix_sub_inplace(int m, int n, double* A, double* B) {
    for (int i = 0; i < m * n; i++) {
        A[i] -= B[i];
    }
}

void matrix_matrix_sub(int m, int n, double* A, double* B, double* C) {
    for (int i = 0; i < m * n; i++) {
        C[i] = A[i] - B[i];
    }
}

void matrix_transpose(int m, int n, double* A, double* B) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            B[i * m + j] = A[j * n + i];
        }
    }
}

void remove_row(int rows, int cols, double* mat, int row_num) {
    for (int i = row_num; i < rows - 1; i++) {
        for (int j = 0; j < cols; j++) {
            mat[i * cols + j] = mat[(i + 1) * cols + j];
        }
    }
}

void remove_cell(int cells, double* vec, int cell_num) {
    for (int i = cell_num; i < cells - 1; i++) {
        vec[i] = vec[i + 1];
    }
}
