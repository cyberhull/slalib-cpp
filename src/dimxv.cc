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
 * Procedure that performs 3D backward unitary transformation (double precision):
 *
 *   vector vb = (inverse of matrix mat) * vector va
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param mat Input matrix; must be unitary, as this routine assumes that the inverse and transpose are identical.
 * @param va Input vector; may be the same as output.
 * @param vb Output: vector; may be the same as input.
 */
void dimxv(const Matrix<double> mat, const Vector<double> va, Vector<double> vb) {
    // inverse of matrix mat * vector va -> vector result
    Vector<double> result;
    for (int j = 0; j < 3; j++) {
        double element = 0.0;
        for (int i = 0; i < 3; i++) {
            element += mat[i][j] * va[i];
        }
        result[j] = element;
    }
    // vector result -> vector vb
    for (int k = 0; k < 3; k++) {
        vb[k] = result[k];
    }
}

}
