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
 * Converts spherical coordinates to direction cosines (double precision).
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param spherical Spherical coordinates in radians: RA/Dec, longitude (positive anticlockwise looking from the
 * positive latitude pole)/latitude, etc.
 * @param cartesian Output: 3-component unit vector; these coordinates are right handed, with the x axis at zero
 *  longitude and latitude, and the z axis at the positive latitude pole.
 */
void dcs2c(const Spherical<double>& spherical, Vector<double> cartesian) {
    const double cos_b = std::cos(spherical.s_b);
    cartesian[0] = std::cos(spherical.s_a) * cos_b;
    cartesian[1] = std::sin(spherical.s_a) * cos_b;
    cartesian[2] = std::sin(spherical.s_b);
}

}
