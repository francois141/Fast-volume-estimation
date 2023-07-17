#ifndef LINALG_H_
#define LINALG_H_
//
// Created by francois on 30.03.23.
//
#include <assert.h>
#include <cblas.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lapacke.h"

/**
 * @brief Wrapper to perform matrix matrix operation, without any transposition
 * or conjugate matrix
 * https://netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gaeda3cbd99c8fb834a60a6412878226e1.html
 *
 * c <- alpha * A * B + beta * C
 *
 * @param[in]  M        Number of rows for the matrix A.
 * @param[in]  N        Number of columns of the matrix B.
 * @param[in]  K        Number of columns of the matrix A.
 * @param[in]  ALPHA    Specifies the scalar alpha
 * @param[in]  BETA     Specifies the scalar beta
 * @param[in]  A        Double precision array. It represents matrix A
 * @param[in]  B        Double precision array. It represents matrix B
 * @param[in, out]  C        Double precision array. It represents matrix C
 *
 * @return void
 */
void matrix_multiplication_wrapper(int m, int n, int k, double alpha,
                                   double beta, double* A, double* B,
                                   double* C);

/**
 * @brief Wrapper to perform the cholesky factorization
 * https://netlib.org/lapack/explore-html/d1/df5/group__real_o_t_h_e_rcomputational_ga2bea5d870e14f7656183aed657af355d.html#ga2bea5d870e14f7656183aed657af355d
 *
 * Finds U, such that A = U * U^T
 *
 * @param[in]      N    Order of the matrix
 * @param[in,out]  A    Input matrix. Returns the decomposition of the matrix
 * (lower triangular)
 *
 * @return 0 on success
 * @return -1 illegal argument
 * @return >0 The matrix is not positive definite
 */
int cholesky_decomposition(int N, double* A);

/**
 * @brief Multiply matrix A with constant c and store the result in A
 *
 * @param[in]      m    Number of rows of the matrix
 * @param[in]      n    Number of columns of the matrix
 * @param[in,out]  A    Input Matrix, size m x n. !Also returns the matrix at
 * the end
 * @param[in]      c    The constant to multiply the matrix with
 */
void matrix_multiply_constant_inplace(int m, int n, double* A, double c);

/**
 * @brief Divide matrix A by constant c and store the result in A
 *
 * @param[in]      m    Number of rows of the matrix
 * @param[in]      n    Number of columns of the matrix
 * @param[in,out]  A    Input Matrix, size m x n. !Also returns the matrix at
 * the end
 * @param[in]      c    The constant to divide the matrix by
 */
void matrix_divide_constant_inplace(int m, int n, double* A, double c);

/**
 * @brief Multiply matrix A with constant c and store the result in matrix B
 *
 * @param[in]      m    Number of rows of the matrix
 * @param[in]      n    Number of columns of the matrix
 * @param[in]      A    Input Matrix, size m x n.
 * @param[in]      c    The constant to multiply the matrix with
 * @param[out]     B    Output Matrix, size m x n. Contains the result at the
 * end
 */
void matrix_multiply_constant(int m, int n, double* A, double c, double* B);

/**
 * @brief Add matrix A with matrix B and store the result in matrix A
 *
 * @param[in]      m    Number of rows of the matrices
 * @param[in]      n    Number of columns of the matrices
 * @param[in,out]  A    Input Matrix, size m x n.
 * @param[in]      B    Input Matrix, size m x n.
 */
void matrix_matrix_add_inplace(int m, int n, double* A, double* B);

/**
 * @brief Add matrix A with matrix B and store the result in matrix C
 *
 * @param[in]      m    Number of rows of the matrices
 * @param[in]      n    Number of columns of the matrices
 * @param[in]      A    Input Matrix, size m x n.
 * @param[in]      B    Input Matrix, size m x n.
 * @param[out]     C    Output Matrix, size m x n.
 */
void matrix_matrix_add(int m, int n, double* A, double* B, double* C);

/**
 * @brief Subtract matrix B from matrix A and store the result in matrix A
 *
 * @param[in]      m    Number of rows of the matrices
 * @param[in]      n    Number of columns of the matrices
 * @param[in,out]  A    Input Matrix, size m x n.
 * @param[in]      B    Input Matrix, size m x n.
 */
void matrix_matrix_sub_inplace(int m, int n, double* A, double* B);

/**
 * @brief Subtract matrix B from matrix A and store the result in matrix C
 *
 * @param[in]      m    Number of rows of the matrices
 * @param[in]      n    Number of columns of the matrices
 * @param[in]      A    Input Matrix, size m x n.
 * @param[in]      B    Input Matrix, size m x n.
 * @param[out]     C    Output Matrix, size m x n.
 */
void matrix_matrix_sub(int m, int n, double* A, double* B, double* C);

/**
 * @brief Computes the transpose of A into B
 *
 * @param[in]      m    Number of rows of the matrix
 * @param[in]      n    Number of columns of the matrix
 * @param[in]      A    Input Matrix, size m x n.
 * @param[out]     A    Output Matrix, size m x n.
 */
void matrix_transpose(int m, int n, double* A, double* B);

/**
 * @brief Removes an entire row from a matrix. Technically wastes
 * memory because it doesn't resize the memory occupied by the
 * matrix
 *
 * @param[in]      m    Number of rows of the matrix
 * @param[in]      n    Number of columns of the matrix
 * @param[in,out]  A    Input Matrix, size m x n.
 */
void remove_row(int rows, int cols, double* mat, int row_num);

/**
 * @brief Removes a cell from a vector. Technically wastes
 * memory because it doesn't resize the memory occupied by the
 * vector
 *
 * @param[in]      m    Number of rows of the matrix
 * @param[in]      n    Number of columns of the matrix
 * @param[in,out]  A    Input Matrix, size m x n.
 */
void remove_cell(int cells, double* vec, int cell_num);

#endif  // LINALG_H_