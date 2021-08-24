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
 * Calculates product of two 3x3 matrices (single precision):
 *
 *   matrix `mc`  =  matrix `ma`  x  matrix `mb`
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param ma First matrix (first multiplicand).
 * @param mb Second matrix (second multiplicand).
 * @param mc Output: product of the two matrices; may be the same matrix as either `ma` or `mb`.
 */
void mxm(const Matrix<float> ma, const Matrix<float> mb, Matrix<float> mc) {
    // multiply into scratch matrix
    Matrix <float> result;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            float element = 0.0f;
            for (int k = 0; k < 3; k++) {
                element += ma[i][k] * mb[k][j];
            }
            result[i][j] = element;
        }
    }

    // return the result
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            mc[i][j] = result[i][j];
        }
    }
}

}
