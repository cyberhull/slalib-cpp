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
 * Calculates position angle of one celestial direction with respect to another (single precision).
 *
 * The procedure sla::bear() performs an equivalent function except that the points are specified in the form
 * of spherical coordinates.
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param v1 Direction cosines of the first point; does *not* have to be a unit vector.
 * @param v2 Direction cosines of the second point; does *not* have to be a unit vector.
 * @return Bearing (position angle), in radians (in the range +/- Pi), of the point v2 with respect to the
 *   point v1; if v2 is a small distance east of v1, returned bearing is about +pi/2; if the two points are
 *   coincident, zero is returned.
 */
float pav(const Vector<float> v1, const Vector<float> v2) {
    // call the double precision version
    Vector<double> dv1, dv2;
    for (int i = 0; i < 3; i++) {
        dv1[i] = v1[i];
        dv2[i] = v2[i];
    }
    return float(dpav(dv1, dv2));
}

}
