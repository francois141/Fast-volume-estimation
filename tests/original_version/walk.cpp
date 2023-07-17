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

#include "vol.h"

#define PI 3.1415926536

// using namespace vol;

double vol::polytope_orig::walk(int k, polytope_orig p) {
    return walk_with_rand_functions(k, &polytope_orig::randd,
                                    &polytope_orig::randi, p);
}

double vol::polytope_orig::walk_with_rand_functions(
    int k, double (polytope_orig::*rand_double)(double),
    int (polytope_orig::*rand_int)(int), polytope_orig& p) {
    double r, max, min, C = 0;
    int dir = (p.*rand_int)(n);
    vec::col_iterator it, itA, end;

    for (it = x.begin(), end = x.end(); it != end; ++it) {
        C += *it * *it;
    }
    C -= x(dir) * x(dir);
    r = sqrt(r2[k + 1] - C);
    max = r - x(dir), min = -r - x(dir);

    // A(x + t v) <= b
    // Av t <= b - Ax
    vec bound = B[dir] - Ai[dir] * x;
    for (it = bound.begin(), end = bound.end(), itA = A.begin_col(dir);
         it != end; ++it, ++itA) {
        if (*itA > 0) {
            if (*it < max) max = *it;
        } else if (*itA < 0)
            if (*it > min) min = *it;
    }

    double t = x(dir) + (p.*rand_double)(max - min) + min;
    x(dir) = t;

    return (C + t * t);
}
