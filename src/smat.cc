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
 *   A is mat non-singular N x N matrix
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
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param n Number of unknowns and size of each matrix dimension.
 * @param mat The `n` x `n` matrix; after the call, contains inverse matrix (if input matrix is singular, contents
 *   after the call is undefined).
 * @param vec The `n`-long vector; after the call, contains solution vector.
 * @param det Return value: determinant; if input matrix is singular, 0.0f is returned.
 * @param ws Workspace (integer array of `n` elements).
 * @return `true` if input matrix is singular, `false` otherwise.
 */
bool smat(int n, float* mat, float* vec, float& det, int* ws) {
    assert(mat && vec && ws);

    // variable-size matrix accessor
    auto element = [mat, n](int i1, int i2) -> float& {
        assert(i1 < n && i2 < n);
        return mat[i1 * n + i2];
    };

    constexpr float EPSILON = 1.0e-20f;
    bool singular = false;
    det = 1.0f;

    int i, k;
    for (k = 0; k < n; k++) {
        float v_max = std::abs(element(k, k));
        int i_max = k;
        if (k != n) {
            for (i = k + 1; i < n; i++) {
                const float t0 = std::abs(element(i, k));
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
                    const float t1 = element(k, j);
                    element(k, j) = element(i_max, j);
                    element(i_max, j) = t1;
                }
                const float t2 = vec[k];
                vec[k] = vec[i_max];
                vec[i_max] = t2;
                det = -det;
            }
            ws[k] = i_max;
            float a_kk = element(k, k);
            det = det * a_kk;
            if (std::abs(det) < EPSILON) {
                singular = -1;
            } else {
                a_kk = 1.0f / a_kk;
                element(k, k) = a_kk;
                for (j = 0; j < n; j++) {
                    if (j != k) {
                        element(k, j) = element(k, j) * a_kk;
                    }
                }
                const float y_k = vec[k] * a_kk;
                vec[k] = y_k;
                for (i = 0; i < n; i++) {
                    const float a_ik = element(i, k);
                    if (i != k) {
                        for (j = 0; j < n; j++) {
                            if (j != k) {
                                element(i, j) = element(i, j) - a_ik * element(k, j);
                            }
                        }
                        vec[i] = vec[i] - a_ik * y_k;
                    }
                }
                for (i = 0; i < n; i++) {
                    if (i != k) {
                        element(i, k) = -element(i, k) * a_kk;
                    }
                }
            }
        }
    }
    if (singular) {
        det = 0.0;
    } else {
        for (k = 0; k < n; k++) {
            const int np1mk = n - 1 - k;
            const int ki = ws[np1mk];
            if (np1mk != ki) {
                for (i = 0; i < n; i++) {
                    const float t3 = element(i, np1mk);
                    element(i, np1mk) = element(i, ki);
                    element(i, ki) = t3;
                }
            }
        }
    }
    return singular;
}

}
