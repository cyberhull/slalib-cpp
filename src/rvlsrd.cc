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
 * Computes velocity component in a given direction due to the Sun's motion with respect to the dynamical Local
 * Standard of Rest (single precision).
 *
 * The Local Standard of Rest used here is the "dynamical" LSR, a point in the vicinity of the Sun which is in a
 * circular orbit around the Galactic centre. The Sun's motion with respect to the dynamical LSR is called the
 * "peculiar" solar motion.
 *
 *  There is another type of LSR, called a "kinematical" LSR. A kinematical LSR is the mean standard of rest of
 *  specified star catalogues or stellar populations, and several slightly different kinematical LSRs are in use.
 *  The Sun's motion with respect to an agreed kinematical LSR is known as the "standard" solar motion. To obtain a
 *  radial velocity correction with respect to an adopted kinematical LSR use the function sla::rvlsrk().
 *
 * Reference:
 *   Delhaye (1965), in "Stars and Stellar Systems", vol 5, p73.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param ra_dec J2000.0 mean RA,Dec (radians).
 * @return Component of "peculiar" solar motion in direction R2000,D2000 (km/s); the result is positive when the Sun
 *   is receding from the given point on the sky.
 */
float rvlsrd(const SphericalDir<float>& ra_dec) {
    /*
     * Peculiar solar motion from Delhaye 1965: in Galactic Cartesian coordinates (+9,+12,+7) km/s. This corresponds
     * to about 16.6 km/s towards Galactic coordinates L2 = 53 deg, B2 = +25 deg, or RA,Dec 17 49 58.7 +28 07 04 J2000.
     *
     * The solar motion is expressed here in the form of a J2000.0 equatorial Cartesian vector:
     *
     *   va[0] = X = -speed * cos(RA) * cos(Dec)
     *   va[1] = Y = -speed * sin(RA) * cos(Dec)
     *   va[2] = Z = -speed * sin(Dec)
     */
    static const vector<float> va = {+0.63823, +14.58542, -7.80116};

    // convert given J2000 RA,Dec to x,y,z
    vector<float> vb;
    cs2c(ra_dec, vb);

    // compute dot product with solar motion vector
    return vdv(va, vb);
}

}
