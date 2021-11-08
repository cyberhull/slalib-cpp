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
 * From the W and V matrices from the SVD factorisation of a matrix (as obtained from the sla::svd() function), obtains
 * the covariance matrix (double precision).
 *
 * Reference:
 *   Numerical Recipes, section 14.3.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param n Number of rows and columns in matrices W and V.
 * @param np First dimension of array containing matrix V.
 * @param nc First dimension of array to receive CVM.
 * @param w `n`x`n` diagonal matrix W (diagonal elements only).
 * @param v `np`x`np` array containing `n`x`n` orthogonal matrix V.
 * @param ws Return value: workspace.
 * @param cvm `nc`x`nc` array to receive covariance matrix.
 */
void svdcov(int n, int np, int nc, const double* w, const double* v, double* ws, double* cvm) {
    assert(n <= np && n <= nc && w && v && ws && cvm);

    auto v_elem = [v, np](int row, int col) -> const double& {
        assert(row < np && col < np);
        return v[row * np + col];
    };
    auto cvm_elem = [cvm, nc](int row, int col) -> double& {
        assert(row < nc && col < nc);
        return cvm[row * nc + col];
    };

    for (int l = 0; l < n; l++) {
        const double s1 = w[l];
        ws[l] = (s1 != 0.0)? 1.0 / (s1 * s1): 0.0;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double s2 = 0.0;
            for (int k = 0; k < n; k++) {
                s2 += v_elem(i, k) * v_elem(j, k) * ws[k];
            }
            cvm_elem(i, j) = s2;
            cvm_elem(j, i) = s2;
        }
    }
}

}
