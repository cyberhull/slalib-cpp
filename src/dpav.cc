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
 * Calculates position angle of one celestial direction with respect to another (double precision).
 *
 * The procedure sla::dbear() performs an equivalent function except that the points are specified in the form
 * of spherical coordinates.
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param va Direction cosines of the first point; does *not* have to be a unit vector.
 * @param vb Direction cosines of the second point; does *not* have to be a unit vector.
 * @return Bearing (position angle), in radians (in the range +/- Pi), of the point vb with respect to the
 *   point va; if vb is a small distance east of va, returned bearing is about +pi/2; if the two points are
 *   coincident, zero is returned.
 */
double dpav(const Vector<double> va, const Vector<double> vb) {
    // unit vector to point 1
    double x1 = va[0];
    double y1 = va[1];
    double z1 = va[2];
    const double length = std::sqrt(x1 * x1 + y1 * y1 + z1 * z1);
    if (length != 0.0) {
        x1 = x1 / length;
        y1 = y1 / length;
        z1 = z1 / length;
    }
    // vector to point 2
    const double x2 = vb[0];
    const double y2 = vb[1];
    const double z2 = vb[2];

    // position angle
    const double sq = y2 * x1 - x2 * y1;
    double cq = z2 * (x1 * x1 + y1 * y1) - z1 * (x2 * x1 + y2 * y1);
    if (sq == 0.0 && cq == 0.0) {
        cq = 1.0;
    }
    return std::atan2(sq, cq);
}

}
