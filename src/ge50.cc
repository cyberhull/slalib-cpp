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
 * Transforms from IAU 1958 galactic coordinates to B1950.0 'FK4' equatorial coordinates (double precision).
 *
 * Reference:
 *   Blaauw et al, Mon.Not.R.Astron.Soc.,121,123 (1960).
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param gal Galactic longitude and latitude L2,B2 (radians).
 * @param fk4 Return value: B1950.0 'FK4' RA,Dec (radians); use the function sla::galeq() if conversion to J2000.0
 *   coordinates is required.
 */
void ge50(const Spherical<double>& gal, Spherical<double>& fk4) {
    /*
     * L2,B2 system of galactic coordinates.
     *
     * P = 192.25    RA of galactic north pole (mean B1950.0; degrees)
     * Q =  62.6     inclination of galactic to mean B1950.0 equator  (degrees)
     * R =  33       longitude of ascending node (degrees)
     *
     * Equatorial to galactic rotation matrix
     *
     * The Euler angles are P, Q, 90-R, about the z then y then z axes.
     *
     *   +CP.CQ.SR-SP.CR    +SP.CQ.SR+CP.CR    -SQ.SR
     *
     *   -CP.CQ.CR-SP.SR    -SP.CQ.CR+CP.SR    +SQ.CR
     *
     *   +CP.SQ             +SP.SQ             +CQ
     */
    static const Matrix<double> mat = {
        {-0.066988739415, -0.872755765852, -0.483538914632},
        {+0.492728466075, -0.450346958020, +0.744584633283},
        {-0.867600811151, -0.188374601723, +0.460199784784},
    };

    // spherical to Cartesian
    Vector<double> v1;
    dcs2c(gal, v1);

    // rotate to mean B1950.0
    Vector<double> v2;
    dimxv(mat, v1, v2);

    // Cartesian to spherical
    Spherical<double> loc;
    dcc2s(v2, loc);

    // introduce E-terms
    addet(loc, 1950.0, fk4);

    // express in conventional ranges
    fk4.set_ra(dranrm(fk4.get_ra()));
    fk4.set_dec(drange(fk4.get_dec()));
}

}
