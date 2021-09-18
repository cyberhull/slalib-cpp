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
 * Transforms from J2000.0 equatorial coordinates to IAU 1958 galactic coordinates (double precision).
 *
 * The equatorial coordinates are J2000.0. Use the function sla::eg50() if conversion from B1950.0 'FK4' coordinates
 * is required.
 *
 * Reference:
 *   Blaauw et al, Mon.Not.R.Astron.Soc.,121,123 (1960).
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param dir J2000.0 RA, Dec (radians).
 * @param gal Galactic longitude and latitude L2, B2 (radians).
 */
void eqgal(const Spherical<double>& dir, Spherical<double>& gal) {
    static const Matrix<double> mat = {
        {-0.054875539726, -0.873437108010, -0.483834985808},
        {+0.494109453312, -0.444829589425, +0.746982251810},
        {-0.867666135858, -0.198076386122, +0.455983795705}
    };
    // spherical to Cartesian
    Vector<double> v1;
    dcs2c(dir, v1);

    // equatorial to galactic
    Vector<double> v2;
    dmxv(mat, v1, v2);

    // Cartesian to spherical
    dcc2s(v2, gal);

    // express in conventional ranges
    gal.set_longitude(dranrm(gal.get_longitude()));
    gal.set_latitude(drange(gal.get_latitude()));
}

}
