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
 * Computes Velocity component in a given direction due to the combination of the rotation of the Galaxy and the
 * motion of the Galaxy relative to the mean motion of the local group (single precision).
 *
 * Reference:
 *   IAU Trans 1976, 168, p201.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param ra_dec J2000.0 mean RA,Dec (radians).
 * @return Component of SOLAR motion in direction `ra_dec` (km/s); the result is positive when the Sun is receding
 *   from the given point on the sky.
 */
float rvlg(const SphericalDir<float>& ra_dec) {
    /*
     * Solar velocity due to Galactic rotation and translation
     *
     * Speed = 300 km/s
     *
     * Apex = L2,B2  90deg, 0deg
     *      = RA,Dec  21 12 01.1  +48 19 47  J2000.0
     *
     * This is expressed in the form of a J2000.0 x,y,z vector:
     *
     *   va[0] = X = -speed * cos(RA) * cos(Dec)
     *   va[1] = Y = -speed * sin(RA) * cos(Dec)
     *   va[2] = Z = -speed * sin(Dec)
     */
    static const vector<float> va = {-148.23284, +133.44888, -224.09467};

    // convert given J2000 RA,Dec to x,y,z
    vector<float> vb;
    cs2c(ra_dec, vb);

    // compute dot product with Solar motion vector
    return vdv(va, vb);
}

}
