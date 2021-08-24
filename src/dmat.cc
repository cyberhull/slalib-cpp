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
 * @param n Number of unknowns (range: [1..3]).
 * @param mat Source matrix; after the call, contains inverse matrix.
 * @param vec Known vector; after the call, contains solution vector.
 * @param det Matrix determinant.
 * @param singular Matrix singularity flag ('true' if the matrix is singular, in which case determinant of 0.0 is
 *   returned, and the contents of `mat` is undefined).
 * @param ws Workspace.
 */
void dmat(int n, matrix<double> mat, Vector<double> vec, double& det, bool& singular, int ws[3]) {
    constexpr double EPSILON = 1.0e-20;

    singular = false;
    det = 1.0;

    int i, j, k;
    for (k = 0; k < n; k++) {
        double amx = std::abs(mat[k][k]);
        int imx = k;
        if (k != n) {
            for (i = k + 1; i < n; i++) {
                const double t0 = std::abs(mat[i][k]);
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
                    const double t1 = mat[k][j];
                    mat[k][j] = mat[imx][j];
                    mat[imx][j] = t1;
                }
                const double t2 = vec[k];
                vec[k] = vec[imx];
                vec[imx] = t2;
                det = -det;
            }
            ws[k] = imx;
            double akk = mat[k][k];
            det = det * akk;
            if (std::abs(det) < EPSILON) {
                singular = true;
            } else {
                akk = 1.0 / akk;
                mat[k][k] = akk;
                for (j = 0; j < n; j++) {
                    if (j != k) {
                        mat[k][j] = mat[k][j] * akk;
                    }
                }
                const double yk = vec[k] * akk;
                vec[k] = yk;
                for (i = 0; i < n; i++) {
                    const double aik = mat[i][k];
                    if (i != k) {
                        for (j = 0; j < n; j++) {
                            if (j != k) {
                                mat[i][j] = mat[i][j] - aik * mat[k][j];
                            }
                        }
                        vec[i] = vec[i] - aik * yk;
                    }
                }
                for (i = 0; i < n; i++) {
                    if (i != k) {
                        mat[i][k] = -mat[i][k] * akk;
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
                    const double t3 = mat[i][np1mk];
                    mat[i][np1mk] = mat[i][ki];
                    mat[i][ki] = t3;
                }
            }
        }
    }
}

}
