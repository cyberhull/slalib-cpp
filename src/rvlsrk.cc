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
 * Calculates velocity component in a given direction due to the Sun's motion with respect to an adopted kinematic
 * Local Standard of Rest (single precision).
 *
 * The Local Standard of Rest used here is one of several "kinematical" LSRs in common use. A kinematical LSR is the
 * mean standard of rest of specified star catalogues or stellar populations. The Sun's motion with respect to a
 * kinematical LSR is known as the "standard" solar motion.
 *
 * There is another sort of LSR, the "dynamical" LSR, which is a point in the vicinity of the Sun which is in a
 * circular orbit around the Galactic centre. The Sun's motion with respect to the dynamical LSR is called the
 * "peculiar" solar motion. To obtain a radial velocity correction with respect to the dynamical LSR use the function
 * sla::rvlsrd().
 *
 * Reference:
 *   Delhaye (1965), in "Stars and Stellar Systems", vol 5, p73.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param pos J2000.0 mean RA,Dec (radians).
 * @return Component of "standard" solar motion in direction R2000,D2000 (km/s); the result is positive when the Sun
 *   is receding from the given point on the sky.
 */
float rvlsrk(const SphericalDir<float>& pos) {
    /*
     * Standard solar motion (from Methods of Experimental Physics, ed Meeks, vol 12, part C, sec 6.1.5.2, p281):
     *
     *   20 km/s towards RA 18h Dec +30d (1900).
     *
     * The solar motion is expressed here in the form of a J2000.0 equatorial Cartesian vector:
     *
     *   va[0] = X = -speed * cos(RA) * cos(Dec)
     *   va[1] = Y = -speed * sin(RA) * cos(Dec)
     *   va[2] = Z = -speed * sin(Dec)
     */
    static const Vector<float> va = {-0.29000, +17.31726, -10.00141};

    // convert given J2000 RA,Dec to x,y,z
    Vector<float> vb;
    cs2c(pos, vb);

    // compute dot product with solar motion vector
    return vdv(va, vb);
}

}
