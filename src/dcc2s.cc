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
 * Converts Cartesian to spherical coordinates (single precision)
 *
 * If `cartesian` is zero-length, zero longitude and latitude and B are returned. At either pole, zero longitude is
 * returned.
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param cartesian Vector representing Cartesian coordinates of a point in space; these coordinates are right handed,
 *  with the x axis at zero longitude and latitude, and the z axis at the positive latitude pole.
 * @param spherical Output: structure containing spherical coordinates (in radians) of the given point; these
 *  coordinates are longitude (positive anticlockwise looking from the positive latitude pole) and latitude.
 */
void dcc2s(const vector<double> cartesian, SphericalDir<double>& spherical) {
    const double x = cartesian[0];
    const double y = cartesian[1];
    const double z = cartesian[2];
    const double r = std::sqrt(x * x + y * y);
    spherical.sd_a = (r == 0.0) ? 0.0 : std::atan2(y, x);
    spherical.sd_b = (z == 0.0) ? 0.0 : std::atan2(z, r);
}

}
