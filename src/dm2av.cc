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
 * Given a rotation matrix, determines corresponding axial vector (double precision).
 *
 * The reference frame rotates clockwise as seen looking along the axial vector from the origin. If rmat is null,
 * so is the result.
 *
 * Original FORTRAN code by P.T.Wallace.
 *
 * @param rmat Rotation matrix; describes a rotation about some arbitrary axis, called the Euler axis.
 * @param axis Output: axial vector; has the same direction as the Euler axis, and its magnitude is the amount of
 *  rotation in radians. The magnitude and direction can be separated by means of the routine sla::vn().
 */
void dm2av(const matrix<double> rmat, vector<double> axis) {
    double x = rmat[2][1] - rmat[1][2];
    double y = rmat[0][2] - rmat[2][0];
    double z = rmat[1][0] - rmat[0][1];
    double s2 = std::sqrt(x * x + y * y + z * z);
    if (s2 != 0.0) {
        double c2 = rmat[0][0] + rmat[1][1] + rmat[2][2] - 1.0;
        double phi = std::atan2(s2 / 2.0, c2 / 2.0);
        double f = phi / s2;
        axis[0] = x * f;
        axis[1] = y * f;
        axis[2] = z * f;
    } else {
        axis[0] = 0.0;
        axis[1] = 0.0;
        axis[2] = 0.0;
    }
}

}
