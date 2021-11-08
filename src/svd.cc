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
#include "f77_utils.h"
#include <cmath>

namespace sla {

/**
 * Singular value decomposition (double precision).
 *
 * This function expresses a given matrix A as the product of three matrices U, W, VT:
 *   A = U x W x VT
 *
 * Where
 *   A   is any M (rows) x N (columns) matrix, where M >= N,
 *   U   is an M x N column-orthogonal matrix,
 *   W   is an N x N diagonal matrix with W[I][I] >= 0,
 *   VT  is the transpose of an N x N orthogonal matrix.
 *
 * Note that M and N, above, are the *logical* dimensions of the matrices and vectors concerned, which can be located
 * in arrays of larger *physical* dimensions, given by MP and NP.
 *
 * References:
 *   The algorithm is an adaptation of the routine SVD in the EISPACK library (Garbow et al 1977, EISPACK Guide
 *   Extension, Springer Verlag), which is a FORTRAN 66 implementation of the Algol routine SVD of Wilkinson &
 *   Reinsch 1971 (Handbook for Automatic Computation, vol 2, ed Bauer et al, Springer Verlag). These references give
 *   full details of the algorithm used here. A good account of the use of SVD in least squares problems is given in
 *   Numerical Recipes (Press et al 1986, Cambridge University Press), which includes another variant of the EISPACK
 *   code.
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param m Number of rows in matrix A.
 * @param n Number of columns in matrix A.
 * @param mp Physical dimension (number of rows) of array containing matrix A.
 * @param np Physical dimension (number of columns) of array containing matrix A.
 * @param a Input: `mp`x`np` array containing `m`x`n` matrix A; output: `mp`x`np` array containing `m`x`n`
 *   column-orthogonal matrix U.
 * @param w Output value: `n`x`n` diagonal matrix W (`n` elements in one-dimensional array of length `np`: diagonal
 *   elements only).
 * @param v Output value: `np`x`np` array containing `n`x`n` orthogonal matrix V; note: it contains matrix V, not the
 *   transpose of matrix V (VT).
 * @param ws Workspace (`np`-long array containing `n` elements).
 * @return 0 = OK, -1 = A of wrong shape, >0 = index of W for which convergence failed; ff returned status is greater
 *   than zero, this need not necessarily be treated as a failure; it means that, due to chance properties of the
 *   matrix A, the QR transformation phase of the routine did not fully converge in a predefined number of
 *   iterations, something that very seldom occurs; when this condition does arise, it is possible that the elements
 *   of the diagonal matrix W have not been correctly found; however, in practice the results are likely to be
 *   trustworthy; applications should report the condition as a warning, but then proceed normally.
 */
int svd(int m, int n, int mp, int np, double* a, double* w, double* v, double* ws) {
    assert(m <= mp && n <= np && a && w && v && ws);

    auto a_elem = [a, m, n, np](int row, int col) -> double& {
        assert(row < m && col < n);
        return a[row * np + col];
    };
    auto v_elem = [v, np](int row, int col) -> double& {
        assert(row < np && col < np);
        return v[row * np + col];
    };
    auto negate = [](double& x) {
        x = -x;
    };

    // maximum number of iterations in QR phase
    constexpr int MAX_ITERATIONS = 30;

    // check that the matrix is the right shape
    int status;
    if (m < n) {
        // no: error status
        status = -1;
    } else {
        // yes: preset the status to OK
        status = 0;

        //
        // Householder reduction to bidiagonal form
        // ----------------------------------------

        double s, f, h, c, x, y, z;
        double g = 0.0;
        double scale = 0.0;
        double an = 0.0;
        int i, j, k, l;
        for (i = 0; i < n; i++) {
            l = i + 1;
            ws[i] = scale * g;
            g = 0.0;
            s = 0.0;
            scale = 0.0;
            if (i < m) {
                for (k = i; k < m; k++) {
                    scale += std::fabs(a_elem(k, i));
                }
                if (scale != 0.0) {
                    for (k = i; k < m; k++) {
                        x = a_elem(k, i) / scale;
                        a_elem(k, i) = x;
                        s += x * x;
                    }
                    f = a_elem(i, i);
                    g = -f_sign(std::sqrt(s), f);
                    h = f * g - s;
                    a_elem(i, i) = f - g;
                    if (i != n - 1) {
                        for (j = l; j < n; j++) {
                            s = 0.0;
                            for (k = i; k < m; k++) {
                                s += a_elem(k, i) * a_elem(k, j);
                            }
                            f = s / h;
                            for (k = i; k < m; k++) {
                                a_elem(k, j) += f * a_elem(k, i);
                            }
                        }
                    }
                    for (k = i; k < m; k++) {
                        a_elem(k, i) *= scale;
                    }
                }
            }
            w[i] = scale * g;
            g = 0.0;
            s = 0.0;
            scale = 0.0;
            if (i < m && i != n - 1) {
                for (k = l; k < n; k++) {
                    scale += std::fabs(a_elem(i, k));
                }
                if (scale != 0.0) {
                    for (k = l; k < n; k++) {
                        x = a_elem(i, k) / scale;
                        a_elem(i, k) = x;
                        s += x * x;
                    }
                    f = a_elem(i, l);
                    g = -f_sign(std::sqrt(s), f);
                    h = f * g - s;
                    a_elem(i, l) = f - g;
                    for (k = l; k < n; k++) {
                        ws[k] = a_elem(i, k) / h;
                    }
                    if (i != m - 1) {
                        for (j = l; j < m; j++) {
                            s = 0.0;
                            for (k = l; k < n; k++) {
                                s += a_elem(j, k) * a_elem(i, k);
                            }
                            for (k = l; k < n; k++) {
                                a_elem(j, k) += s * ws[k];
                            }
                        }
                    }
                    for (k = l; k < n; k++) {
                        a_elem(i, k) *= scale;
                    }
                }
            }
            // overestimate of largest column norm for convergence test
            an = std::fmax(an, std::fabs(w[i]) + std::fabs(ws[i]));
        }

        //
        // Accumulation of right-hand transformations
        // ------------------------------------------

        for (i = n - 1; i >= 0; i--) {
            if (i != n - 1) {
                if (g != 0.0) {
                    for (j = l; j < n; j++) {
                        v_elem(j, i) = (a_elem(i, j) / a_elem(i, l)) / g;
                    }
                    for (j = l; j < n; j++) {
                        s = 0.0;
                        for (k = l; k < n; k++) {
                            s += a_elem(i, k) * v_elem(k, j);
                        }
                        for (k = l; k < n; k++) {
                            v_elem(k, j) += s * v_elem(k, i);
                        }
                    }
                }
                for (j = l; j < n; j++) {
                    v_elem(i, j) = 0.0;
                    v_elem(j, i) = 0.0;
                }
            }
            v_elem(i, i) = 1.0;
            g = ws[i];
            l = i;
        }

        //
        // Accumulation of left-hand transformations
        // -----------------------------------------

        for (i = n - 1; i >= 0; i--) {
            l = i + 1;
            g = w[i];
            if (i != n - 1) {
                for (j = l; j < n; j++) {
                    a_elem(i, j) = 0.0;
                }
            }
            if (g != 0.0) {
                if (i != n - 1) {
                    for (j = l; j < n; j++) {
                        s = 0.0;
                        for (k = l; k < m; k++) {
                            s += a_elem(k, i) * a_elem(k, j);
                        }
                        f = (s / a_elem(i, i)) / g;
                        for (k = i; k < m; k++) {
                            a_elem(k, j) += f * a_elem(k, i);
                        }
                    }
                }
                for (j = i; j < m; j++) {
                    a_elem(j, i) /= g;
                }
            } else {
                for (j = i; j < m; j++) {
                    a_elem(j, i) = 0.0;
                }
            }
            a_elem(i, i) += 1.0;
        }

        //
        // Diagonalisation of the bidiagonal form
        // --------------------------------------

        for (k = n - 1; k >= 0; k--) {
            const int k1 = k - 1;

            // iterate until converged
            int iterations = 0;
            while (iterations < MAX_ITERATIONS) {
                iterations++;

                // test for splitting into submatrices
                bool cancel = true;
                int l1;
                for (l = k; l >= 0; l--) {
                    l1 = l - 1;
                    if (an + std::fabs(ws[l]) == an) {
                        cancel = false;
                        break;
                    }
                    // (Following never attempted for l=0 because ws[0] is zero)
                    if (an + std::fabs(w[l1]) == an) {
                        break;
                    }
                }

                // Cancellation of ws[l] if l>1
                if (cancel) {
                    s = 1.0;
                    for (i = l; i <= k; i++) {
                        f = s * ws[i];
                        if (an + std::fabs(f) == an) {
                            break;
                        }
                        g = w[i];
                        h = std::sqrt(f * f + g * g);
                        w[i] = h;
                        c = g / h;
                        s = -f / h;
                        for (j = 0; j < m; j++) {
                            y = a_elem(j, l1);
                            z = a_elem(j, i);
                            a_elem(j, l1) = y * c + z * s;
                            a_elem(j, i) = -y * s + z * c;
                        }
                    }
                }

                // converged?
                z = w[k];
                if (l == k) {

                    // yes:  stop iterating
                    iterations = MAX_ITERATIONS;

                    // ensure singular values are non-negative
                    if (z < 0.0) {
                        w[k] = -z;
                        for (j = 0; j < n; j++) {
                            negate(v_elem(j, k));
                        }
                    }
                } else {

                    // not converged yet: set status if iteration limit reached
                    if (iterations == MAX_ITERATIONS) {
                        status = k;
                    }

                    // shift from bottom 2x2 minor
                    x = w[l];
                    y = w[k1];
                    g = ws[k1];
                    h = ws[k];
                    f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
                    const double abs_f = std::fabs(f);
                    g = abs_f > 1.0e15? abs_f: std::sqrt(f * f + 1.0);
                    f = ((x - z) * (x + z) + h * (y / (f + f_sign(g, f)) - h)) / x;

                    // next QR transformation
                    c = 1.0;
                    s = 1.0;
                    for (int i1 = l; i1 <= k1; i1++) {
                        i = i1 + 1;
                        g = ws[i];
                        y = w[i];
                        h = s * g;
                        g = c * g;
                        z = std::sqrt(f * f + h * h);
                        ws[i1] = z;
                        if (z != 0.0) {
                            c = f / z;
                            s = h / z;
                        } else {
                            c = 1.0;
                            s = 0.0;
                        }
                        f = x * c + g * s;
                        g = -x * s + g * c;
                        h = y * s;
                        y = y * c;
                        for (j = 0; j < n; j++) {
                            x = v_elem(j, i1);
                            z = v_elem(j, i);
                            v_elem(j, i1) = x * c + z * s;
                            v_elem(j, i) = -x * s + z * c;
                        }
                        z = std::sqrt(f * f + h * h);
                        w[i1] = z;
                        if (z != 0.0) {
                            c = f / z;
                            s = h / z;
                        }
                        f = c * g + s * y;
                        x = -s * g + c * y;
                        for (j = 0; j < m; j++) {
                            y = a_elem(j, i1);
                            z = a_elem(j, i);
                            a_elem(j, i1) = y * c + z * s;
                            a_elem(j, i) = -y * s + z * c;
                        }
                    }
                    ws[l] = 0.0;
                    ws[k] = f;
                    w[k] = x;
                }
            }
        }
    }
    return status;
}

}
