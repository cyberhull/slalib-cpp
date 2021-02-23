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
 * Procedure that forms rotation matrix corresponding to a given axial vector (single precision).
 *
 * The reference frame rotates clockwise as seen looking along the axial vector from the origin.
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param vec Axial vector (in radians); has the same direction as the Euler axis, and its magnitude is the
 *  amount of rotation in radians.
 * @param mat Output: rotation matrix that describes a rotation about some arbitrary axis, called the Euler axis. If
 *   `vec` is zero-length, the unit matrix is returned.
 */
void av2m(const vector<float> vec, matrix<float> mat) {
    // rotation angle and magnitude of axial vector
    float x = vec[0];
    float y = vec[1];
    float z = vec[2];
    const float phi = std::sqrt(x * x + y * y + z * z);
    const float sin_phi = std::sin(phi);
    const float cos_phi = std::cos(phi);
    const float w = 1.0f - cos_phi;

    // Euler axis - direction of axial vector
    if (phi != 0.0f) {
        x = x / phi;
        y = y / phi;
        z = z / phi;
    }
    // compute and return the rotation matrix
    mat[0][0] = x * x * w + cos_phi;
    mat[0][1] = x * y * w + z * sin_phi;
    mat[0][2] = x * z * w - y * sin_phi;
    mat[1][0] = x * y * w - z * sin_phi;
    mat[1][1] = y * y * w + cos_phi;
    mat[1][2] = y * z * w + x * sin_phi;
    mat[2][0] = x * z * w + y * sin_phi;
    mat[2][1] = y * z * w - x * sin_phi;
    mat[2][2] = z * z * w + cos_phi;
}

}
