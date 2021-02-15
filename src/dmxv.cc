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
 * Procedure that performs the 3D forward unitary transformation (double precision):
 *
 *   vector `vb` = matrix `rm` * vector `va`
 *
 * Original FORTRAN code by P.T.Wallace.
 *
 * @param rm Transformation matrix.
 * @param va Vector to transform.
 * @param vb Output: vector `va` multiplied by matrix `rm`; can be the same as `va`.
 */
void dmxv(const matrix<double> rm, const vector<double> va, vector<double> vb) {
    // matrix `rm` * vector `va` -> vector `result`
    vector<float> result;
    for (int j = 0; j < 3; j++) {
        double element = 0.0f;
        for (int i = 0; i < 3; i++) {
            element += rm[i][j] * va[i];
        }
        result[j] = element;
    }

    // return result vector
    for (int i = 0; i < 3; i++) {
        vb[i] = result[i];
    }
}

}
