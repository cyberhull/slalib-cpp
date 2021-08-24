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
 * Transforms from de Vaucouleurs supergalactic coordinates to IAU 1958 galactic coordinates (double precision).
 *
 * References (these two references give different values for the galactic longitude of the supergalactic origin;
 * both are wrong; the correct value is L2=137.37.):
 *
 *   de Vaucouleurs, de Vaucouleurs, & Corwin, Second Reference Catalogue of Bright Galaxies, U. Texas, page 8.
 *
 *   Systems & Applied Sciences Corp., Documentation for the machine-readable version of the above catalogue,
 *   Contract NAS 5-26490.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param sgal Supergalactic longitude and latitude (radians).
 * @param gal Return value: galactic longitude and latitude L2,B2 (radians).
 */
void supgal(const Spherical<double>& sgal, Spherical<double>& gal) {
    /*
     *  System of supergalactic coordinates:
     *
     *    SGL   SGB        L2     B2      (deg)
     *     -    +90      47.37  +6.32
     *     0     0         -      0
     *
     * Galactic to supergalactic rotation matrix:
     */
    static const Matrix<double> mat = {
        {-0.735742574804, +0.677261296414, +0.000000000000},
        {-0.074553778365, -0.080991471307, +0.993922590400},
        {+0.673145302109, +0.731271165817, +0.110081262225}
    };
    // spherical to Cartesian
    Vector<double> v1;
    dcs2c(sgal, v1);

    // supergalactic to galactic
    Vector<double> v2;
    dimxv(mat, v1, v2);

    // Cartesian to spherical
    dcc2s(v2, gal);

    // express in conventional ranges
    gal.set_longitude(dranrm(gal.get_longitude()));
    gal.set_latitude(drange(gal.get_latitude()));
}

}
