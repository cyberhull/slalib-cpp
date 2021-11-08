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

namespace sla {

/**
 * From a given vector and the SVD of a matrix (as obtained from the sla::svd() function), obtains the solution
 * vector (double precision).
 *
 * This function solves the equation:
 *   A . x = b
 *
 * where
 *   A  is a given M (rows) x N (columns) matrix, where M >= N,
 *   x  is the N-vector we wish to find,
 *   b  is a given M-vector
 *
 * by means of the Singular Value Decomposition method (SVD). In this method, the matrix A is first factorised (for
 * example by the function sla::svd()) into the following components:
 *   A = U x W x VT
 *
 * where:
 *   A   is the M (rows) x N (columns) matrix,
 *   U   is an M x N column-orthogonal matrix,
 *   W   is an N x N diagonal matrix with W[I][I] >= 0
 *   VT  is the transpose of an NxN orthogonal matrix.
 *
 * Note that M and N, above, are the *logical* dimensions of the matrices and vectors concerned, which can be
 * located in arrays of larger *physical* dimensions MP and NP.
 *
 * The solution is found from the expression:
 *    x = V . [diag(1/Wj)] . (transpose(U) . b)
 *
 * If matrix A is square, and if the diagonal matrix W is not adjusted, the method is equivalent to conventional
 * solution of simultaneous equations. If M>N, the result is a least-squares fit.
 *
 * If the solution is poorly determined, this shows up in the SVD factorisation as very small or zero Wj values. Where
 * a Wj value is small but non-zero it can be set to zero to avoid ill effects. The present function detects such zero
 * Wj values and produces a sensible solution, with highly correlated terms kept under control rather than being
 * allowed to elope to infinity, and with meaningful values for the other terms.
 *
 * Reference:
 *   Numerical Recipes, section 2.9.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param m Numbers of rows in matrix A.
 * @param n Numbers of columns in matrix A.
 * @param mp Physical dimension (number of rows) of array containing matrix A.
 * @param np Physical dimension (number of columns) of array containing matrix A.
 * @param b Known vector.
 * @param u `mp`x`np` array containing `m`x`n` matrix U.
 * @param w `n`x`n` diagonal matrix W (diagonal elements only).
 * @param v `np`x`np` array containing `n`x`n` orthogonal matrix V.
 * @param ws Return value: workspace (`n`-element vector).
 * @param x Return value: unknown `n`-element vector x.
 */
void svdsol(int m, int n, int mp, int np, const double* b, const double* u, const double* w, const double* v,
    double* ws, double* x) {
    assert(m <= mp && n <= np && b && u && w && v && ws && x);

    auto u_elem = [u, m, n, np](int row, int col) -> const double& {
        assert(row < m && col < n);
        return u[row * np + col];
    };
    auto v_elem = [v, np](int row, int col) -> const double& {
        assert(row < np && col < np);
        return v[row * np + col];
    };

    // calculate [diag(1/Wj)] . transpose(u) . b (or zero for zero Wj)
    for (int j = 0; j < n; j++) {
        double s1 = 0.0;
        if (w[j] != 0.0) {
            for (int i = 0; i < m; i++) {
                s1 += u_elem(i, j) * b[i];
            }
            s1 = s1 / w[j];
        }
        ws[j] = s1;
    }

    // multiply by matrix v to get result
    for (int k = 0; k < n; k++) {
        double s2 = 0.0;
        for (int l = 0; l < n; l++) {
            s2 += v_elem(k, l) * ws[l];
        }
        x[k] = s2;
    }
}

}
