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
 * Calculates angle between two points on sa sphere (double precision).
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param sa Spherical coordinates of one point (RA, longitude, etc.) (radians).
 * @param sb Spherical coordinates of the other point (DEC, latitude, etc.) (radians).
 * @return The angle, in radians, between the two points; it is always positive.
 */
double dsep(const Spherical<double>& sa, const Spherical<double>& sb) {
    // convert coordinates from spherical to Cartesian
    Vector<double> va, vb;
    dcs2c(sa, va);
    dcs2c(sb, vb);

    // calculate and return angle between the vectors
    return dsepv(va, vb);
}

}
