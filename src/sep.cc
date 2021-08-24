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
 * Calculates angle between two points on sa sphere (single precision).
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param sa Spherical coordinates of one point (RA, longitude, etc.) (radians).
 * @param sb Spherical coordinates of the other point (DEC, latitude, etc.) (radians).
 * @return The angle, in radians, between the two points; it is always positive.
 */
float sep(const SphericalDir<float>& sa, const SphericalDir<float>& sb) {
    // use double precision version
    const SphericalDir<double> da = {sa.sd_a, sa.sd_b};
    const SphericalDir<double> db = {sb.sd_a, sb.sd_b};
    return (float) dsep(da, db);
}

}
