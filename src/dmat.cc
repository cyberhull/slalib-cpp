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
 * Matrix inversion and solution of simultaneous equations (double precision).
 *
 * For the set of n simultaneous equations in n unknowns:
 *    M . V = X
 *
 * where:
 *    M is a non-singular N x N matrix
 *    V is the vector of N unknowns
 *    X is the known vector
 *
 * sla::dmat() computes:
 *    the inverse of matrix M
 *    the determinant of matrix M
 *    the vector of N unknowns
 *
 * Algorithm is Gaussian elimination with partial pivoting. Speed: very fast. The function is fairly accurate - errors
 * 1 to 4 times those of functions optimized for accuracy.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param n Number of unknowns.
 * @param mat Source matrix of size [n][n]; after the call, contains inverse matrix (if source matrix is singular,
 *   contents after the call is undefined).
 * @param vec Known vector of `n` elements; after the call, contains solution vector.
 * @param det Return value: determinant; if input matrix is singular, 0.0 is returned.
 * @param ws Workspace (integer array of `n` elements).
 * @return `true` if input matrix is singular, `false` otherwise.
 */
bool dmat(int n, double* mat, double* vec, double& det, int* ws) {
    assert(mat && vec && ws);

    // variable-size matrix accessor
    auto element = [mat, n](int i1, int i2) -> double& {
        assert(i1 < n && i2 < n);
        return mat[i1 * n + i2];
    };

    constexpr double EPSILON = 1.0e-20;
    bool singular = false;
    det = 1.0;

    int i, j, k;
    for (k = 0; k < n; k++) {
        double amx = std::abs(element(k, k));
        int imx = k;
        if (k != n) {
            for (i = k + 1; i < n; i++) {
                const double t0 = std::abs(element(i, k));
                if (t0 > amx) {
                    amx = t0;
                    imx = i;
                }
            }
        }
        if (amx < EPSILON) {
            singular = true;
        } else {
            if (imx != k) {
                for (j = 0; j < n; j++) {
                    const double t1 = element(k, j);
                    element(k, j) = element(imx, j);
                    element(imx, j) = t1;
                }
                const double t2 = vec[k];
                vec[k] = vec[imx];
                vec[imx] = t2;
                det = -det;
            }
            ws[k] = imx;
            double akk = element(k, k);
            det = det * akk;
            if (std::abs(det) < EPSILON) {
                singular = true;
            } else {
                akk = 1.0 / akk;
                element(k, k) = akk;
                for (j = 0; j < n; j++) {
                    if (j != k) {
                        element(k, j) = element(k, j) * akk;
                    }
                }
                const double yk = vec[k] * akk;
                vec[k] = yk;
                for (i = 0; i < n; i++) {
                    const double aik = element(i, k);
                    if (i != k) {
                        for (j = 0; j < n; j++) {
                            if (j != k) {
                                element(i, j) = element(i, j) - aik * element(k, j);
                            }
                        }
                        vec[i] = vec[i] - aik * yk;
                    }
                }
                for (i = 0; i < n; i++) {
                    if (i != k) {
                        element(i, k) = -element(i, k) * akk;
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
                    const double t3 = element(i, np1mk);
                    element(i, np1mk) = element(i, ki);
                    element(i, ki) = t3;
                }
            }
        }
    }
    return singular;
}

}
