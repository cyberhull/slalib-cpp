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
 * Calculates bearing (position angle) of one point on a sphere relative to another given spherical coordinates of
 * the point (single precision).
 *
 * Original FORTRAN code by Rutherford Appleton Laboratory / P.T. Wallace.
 *
 * @param a1 RA/longitude/etc. of the first point, radians.
 * @param b1 Dec/latitude/etc. of the first point, radians.
 * @param a2 RA/longitude/etc. of the second point, radians.
 * @param b2 Dec/latitude/etc. of the second point, radians.
 * @return Bearing (position angle), in radians (range range +/- pi), of point a2,b2 as seen from point a1,b1; if
 *  a2,b2 is due east of a1,b1 the bearing is +pi/2; if the two points are coincident, zero is returned.
 */
double bear(float a1, float b1, float a2, float b2) {
    const float cos_b2 = std::cos(b2);
    const float da = a2 - a1;
    const float x = std::sin(b2) * std::cos(b1) - cos_b2 * std::sin(b1) * std::cos(da);
    const float y = std::sin(da) * cos_b2;
    return (x != 0.0f || y != 0.0f) ? std::atan2(y, x) : 0.0f;
}

}
