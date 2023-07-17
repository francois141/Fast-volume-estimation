/***********************************************************************
 *  This code is part of PolyVest.
 *
 *  Copyright (C) 2013, 2016 Cunjing Ge, Institute of Software, Chinese
 *  Academy of Sciences, Beijing, China. All rights reserved.
 *  E-mail: <gecj@ios.ac.cn>.
 *
 *  PolyVest is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PolyVest is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PolyVest. If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/

#include <algorithm>
#include <cmath>
#include <iostream>

#include "armadillo"
#include "memory.h"
#include "time.h"

#ifndef POLYVOL_H
#define POLYVOL_H

using namespace std;
using namespace arma;

namespace vol {

class polytope_orig {
   public:
    polytope_orig(int rows, int cols);
    ~polytope_orig();
    double matA(double val, int i, int j) {
        A(i, j) = val;
        return val;
    }
    double matA(int i, int j) { return A(i, j); }
    double vecb(double val, int i) {
        b(i) = val;
        return val;
    }
    double vecb(int i) { return b(i); }

    bool msg_off;
    bool check_planes_off;
    double walk(int k);
    double walk(int k, polytope_orig p);
    double walk_with_rand_functions(
        int k, double (polytope_orig::*rand_double)(double),
        int (polytope_orig::*rand_int)(int), polytope_orig &p);
    void Preprocess();
    double EstimateVol(int coef);
    double EstimateVol(int coef, double (polytope_orig::*rand_double)(double),
                       int (polytope_orig::*rand_int)(int));
    double Volume() const { return vol; }
    void checkHPs();
    void genInitE(double &R2, vec &Ori);
    double randd(double u) { return rand() * u / RAND_MAX; }
    int randi(int u) { return rand() % u; }

    double randd_test(double u) { return 1.0; }
    int randi_test(int u) { return 1.0; }

    // polytope_orig denoted by: Ax<=b, A is an (m x n) matrix.
    int m, n;
    mat A;
    vec b;

    // approximating volume variables
    vec x;
    double vol, determinant;
    int l;
    double *r2;

    // reciprocal of beta, beta-cut
    double beta_r;

    vec *B;
    mat *Ai;
};

inline polytope_orig::polytope_orig(int rows, int cols)
    : msg_off(false),
      check_planes_off(false),
      m(rows),
      n(cols),
      A(rows, cols, fill::zeros),
      b(rows, fill::zeros),
      x(cols, fill::zeros),
      vol(0),
      determinant(0) {
    srand((int)time(0));

    beta_r = 2 * n;  // 2 * n;

    l = (int)(n * log((double)beta_r) / log((double)2)) + 2;
    r2 = new double[l];
    for (int i = 0; i < l; i++) r2[i] = pow((double)2, (double)(2 * i) / n);

    B = new vec[n];
    Ai = new mat[n];
    rowvec exp(n);
    exp.ones();
    for (int i = 0; i < n; i++) {
        B[i] = b / A.col(i);
        Ai[i] = A / (A.col(i) * exp);
        // std::cout << "size Ai: \n" << arma::size(Ai[i]);
        // std::cout << "size A: \n" << arma::size(A);
        // std::cout << "size A col: \n" << arma::size(A.col(i));
    }
}

inline polytope_orig::~polytope_orig() {
    delete[] r2;
    delete[] B;
    delete[] Ai;
}

}  // namespace vol

#endif
