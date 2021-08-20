/*
 * C++ Port of the SLALIB library.
 * Written by Vadim Sytnikov.
 * Copyright (C) 2021 CyberHULL, Ltd.
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 */
#include "slalib.h"
#include <cmath>

namespace sla {

/**
 * Performs matrix inversion & solution of simultaneous equations (single precision).
 *
 * For the set of n simultaneous equations in n unknowns:
 *   A . Y = X
 *
 * where:
 *   A is a non-singular N x N matrix
 *   Y is the vector of N unknowns
 *   X is the known vector
 *
 * sla::smatrx() computes:
 *    the inverse of matrix A
 *    the determinant of matrix A
 *    the vector of N unknowns
 *
 * Algorithm: Gaussian elimination with partial pivoting. Speed: very fast. Accuracy: fairly accurate - errors 1 to 4
 * times those of routines optimised for accuracy.
 *
 * Original FORTRAN code by Rutherford Appleton Laboratory / P.T. Wallace.
 *
 * @param n Number of unknowns and size of each matrix dimension.
 * @param a The `n` x `n` matrix; after the call, contains inverse matrix (if input matrix is singular, contents
 *   after the call is undefined).
 * @param y The `n`-long vector; after the call, contains solution vector.
 * @param d Return value: determinant; if input matrix is singular, 0.0f is returned.
 * @param w Workspace.
 * @return `true` if input matrix is singular, `false` otherwise.
 */
bool smat(int n, float* a, float* y, float& d, int* w) {
    constexpr float EPSILON = 1.0e-20f;

    // variable-size matrix accessor
    auto mat = [a, n](int i1, int i2) -> float& {
        return a[i2 * n + i1];
    };

    bool singular = false;
    d = 1.0f;
    int i, k;
    for (k = 0; k < n; k++) {
        float v_max = std::abs(mat(k, k));
        int i_max = k;
        if (k != n) {
            for (i = k + 1; i < n; i++) {
                const float t0 = std::abs(mat(i, k));
                if (t0 > v_max) {
                    v_max = t0;
                    i_max = i;
                }
            }
        }
        if (v_max < EPSILON) {
            singular = true;
        } else {
            int j;
            if (i_max != k) {
                for (j = 0; j < n; j++) {
                    const float t1 = mat(k, j);
                    mat(k, j) = mat(i_max, j);
                    mat(i_max, j) = t1;
                }
                const float t2 = y[k];
                y[k] = y[i_max];
                y[i_max] = t2;
                d = -d;
            }
            w[k] = i_max;
            float a_kk = mat(k, k);
            d = d * a_kk;
            if (std::abs(d) < EPSILON) {
                singular = -1;
            } else {
                a_kk = 1.0f / a_kk;
                mat(k, k) = a_kk;
                for (j = 0; j < n; j++) {
                    if (j != k) {
                        mat(k, j) = mat(k, j) * a_kk;
                    }
                }
                const float y_k = y[k] * a_kk;
                y[k] = y_k;
                for (i = 0; i < n; i++) {
                    const float a_ik = mat(i, k);
                    if (i != k) {
                        for (j = 0; j < n; j++) {
                            if (j != k) {
                                mat(i, j) = mat(i, j) - a_ik * mat(k, j);
                            }
                        }
                        y[i] = y[i] - a_ik * y_k;
                    }
                }
                for (i = 0; i < n; i++) {
                    if (i != k) {
                        mat(i, k) = -mat(i, k) * a_kk;
                    }
                }
            }
        }
    }
    if (singular) {
        d = 0.0;
    } else {
        for (k = 0; k < n; k++) {
            const int np1mk = n - 1 - k;
            const int ki = w[np1mk];
            if (np1mk != ki) {
                for (i = 0; i < n; i++) {
                    const float t3 = mat(i, np1mk);
                    mat(i, np1mk) = mat(i, ki);
                    mat(i, ki) = t3;
                }
            }
        }
    }
    return singular;
}

}
