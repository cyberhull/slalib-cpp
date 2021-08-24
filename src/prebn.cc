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
 * Generates the matrix of precession between two epochs, using the old, pre-IAU1976, Bessel-Newcomb model, using
 * Kinoshita's formulation (double precision).
 *
 * Reference:
 *   Kinoshita, H. (1975) 'Formulas for precession', SAO Special Report No. 364, Smithsonian Institution Astrophysical
 *   Observatory, Cambridge, Massachusetts.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param be0 Beginning Besselian epoch.
 * @param be1 Ending Besselian epoch.
 * @param mat Return value: precession matrix; the matrix is in the sense v(be1) = mat * v(ep0).
 */
void prebn(double be0, double be1, Matrix<double> mat) {
    // arc seconds to radians
    constexpr double ARCSECS_2_RADIANS = 0.484813681109535994e-5;

    // interval between basic epoch B1850.0 and beginning epoch in TC
    const double bigt = (be0 - 1850.0) / 100.0;

    // interval over which precession required, in tropical centuries
    const double t = (be1 - be0) / 100.0;

    // Euler angles
    const double tas2r = t * ARCSECS_2_RADIANS;
    const double w = 2303.5548 + (1.39720 + 0.000059 * bigt) * bigt;

    const double zeta = (w + (0.30242 - 0.000269 * bigt + 0.017996 * t) * t) * tas2r;
    const double z = (w + (1.09478 + 0.000387 * bigt + 0.018324 * t) * t) * tas2r;
    const double theta = (2005.1125 + (-0.85294 - 0.000365 * bigt) * bigt +
        (-0.42647 - 0.000365 * bigt - 0.041802 * t) * t) * tas2r;

    // rotation matrix
    deuler("ZYZ", -zeta, theta, -z, mat);
}

}
