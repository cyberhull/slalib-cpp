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
 * Calculates velocity component in a given direction due to the rotation of the Galaxy (single precision).
 *
 * The Local Standard of Rest used here is a point in the vicinity of the Sun which is in a circular orbit around the
 * Galactic centre. Sometimes called the "dynamical" LSR, it is not to be confused with a "kinematical" LSR, which is
 * the mean standard of rest of star catalogues or stellar populations.
 *
 * Reference: The orbital speed of 220 km/s used here comes from Kerr & Lynden-Bell (1986), MNRAS, 221, p1023.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param pos J2000.0 mean RA,Dec (radians).
 * @return Component of dynamical LSR motion in direction R2000,D2000 (km/s); the result is positive when the
 *   dynamical LSR is receding from the given point on the sky.
 */
float rvgalc(const SphericalDir<float>& pos) {
    /*
     * LSR velocity due to Galactic rotation
     *
     * Speed = 220 km/s
     * Apex  = L2,B2  90deg, 0deg
     *       = RA,Dec  21 12 01.1  +48 19 47  J2000.0
     *
     * This is expressed in the form of a J2000.0 x,y,z vector:
     *   va[0] = X = -speed * cos(RA) * cos(Dec)
     *   va[1] = Y = -speed * sin(RA) * cos(Dec)
     *   va[2] = Z = -speed * sin(Dec)
     */
    static const Vector<float> va = {-108.70408, +97.86251, -164.33610};

    // convert given J2000 RA,Dec to x,y,z
    Vector<float> vb;
    cs2c(pos, vb);

    // compute dot product with LSR motion vector
    return vdv(va, vb);
}

}
