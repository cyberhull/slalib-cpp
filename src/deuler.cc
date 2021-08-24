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
#include <cstring>

namespace sla {

/**
 * Form a rotation matrix from the Euler angles - three successive rotations about specified
 * Cartesian axes (double precision).
 *
 * A rotation is positive when the reference frame rotates anticlockwise as seen looking towards the origin from the
 * positive region of the specified axis.
 *
 * The characters of `order` define which axes the three successive rotations are about. A typical value is 'ZXZ',
 * indicating that `mat` is to become the direction cosine matrix corresponding to rotations of the reference frame
 * through `phi` radians about the old Z-axis, followed by `theta` radians about the resulting X-axis, then `psi`
 * radians about the resulting Z-axis.
 *
 * The axis names can be any of the following, in any order or combination: X, Y, Z, uppercase or lowercase,
 * 1, 2, 3. Normal axis labelling/numbering conventions apply; the xyz (=123) triad is right-handed. Thus, the 'ZXZ'
 * example given above could be written 'zxz' or '313' (or even 'ZxZ' or '3xZ'). `order` is terminated by length or
 * by the first unrecognized character.
 *
 * Fewer than three rotations are acceptable, in which case the later angle arguments are ignored. If all rotations
 * are zero, the identity matrix is produced.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param order Specifies about which axes the rotations occur.
 * @param phi First rotation (radians).
 * @param theta Second rotation (radians).
 * @param psi Third rotation (radians).
 * @param mat Output: rotation matrix.
 */
void deuler(const char* order, const double phi, const double theta, const double psi, Matrix<double> mat) {
    // set result matrix to identity matrix
    Matrix<double> result;
    int i, j;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            result[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // establish length of axis string
    int l = std::strlen(order);
    if (l > 3) {
        l = 3;
    }

    // look at each character of axis string until finished
    for (int n = 0; n < l; n++) {

        // initialize rotation matrix for the current rotation
        Matrix<double> rotation;
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                rotation[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }

        // pick up the appropriate Euler angle and take sine & cosine
        double angle;
        switch (n) {
            case 0:
                angle = phi;
                break;
            case 1:
                angle = theta;
                break;
            default:
                angle = psi;
        }
        const double sin_a = sin(angle);
        const double cos_a = cos(angle);

        // identify the axis
        const char axis = order[n];
        switch (axis) {
            case 'X':
            case 'x':
            case '1':
                // matrix for x-rotation
                rotation[1][1] = cos_a;
                rotation[1][2] = sin_a;
                rotation[2][1] = -sin_a;
                rotation[2][2] = cos_a;
                break;
            case 'Y':
            case 'y':
            case '2':
                // matrix for y-rotation
                rotation[0][0] = cos_a;
                rotation[0][2] = -sin_a;
                rotation[2][0] = sin_a;
                rotation[2][2] = cos_a;
                break;
            case 'Z':
            case 'z':
            case '3':
                // matrix for z-rotation
                rotation[0][0] = cos_a;
                rotation[0][1] = sin_a;
                rotation[1][0] = -sin_a;
                rotation[1][1] = cos_a;
                break;
            default:
                // unrecognized character - fake end of string
                l = 0;
        }

        // apply the current rotation (matrix rotation x matrix result)
        Matrix<double> combined;
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                double element = 0.0;
                for (int k = 0; k < 3; k++) {
                    element += rotation[i][k] * result[k][j];
                }
                combined[i][j] = element;
            }
        }
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                result[i][j] = combined[i][j];
            }
        }
    }

    // copy the result
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            mat[i][j] = result[i][j];
        }
    }
}

}
