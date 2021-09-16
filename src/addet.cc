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
 * Adds the E-terms (elliptic component of annual aberration) to a pre IAU 1976 mean place to conform to the old
 * catalogue convention (double precision).
 *
 * Most star positions from pre-1984 optical catalogues (or derived from astrometry using such stars) embody the
 * E-terms. If it is necessary to convert a formal mean place (for example a pulsar timing position) to one consistent
 * with such a star catalogue, then the RA,Dec should be adjusted using this routine.
 *
 * Reference:
 *   Explanatory Supplement to the Astronomical Ephemeris, section 2D, page 48.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param dir RA,Dec without E-terms (radians).
 * @param be Besselian epoch of mean equator and equinox.
 * @param edir Return value: RA,Dec with E-terms included (radians).
 */
void addet(const Spherical<double>& dir, double be, Spherical<double>& edir) {
    // get E-terms vector
    Vector<double> et;
    etrms(be, et);

    // spherical to Cartesian
    Vector<double> v;
    dcs2c(dir, v);

    // include the E-terms
    for (int i = 0; i < 3; i++) {
        v[i] += et[i];
    }
    // Cartesian to spherical
    dcc2s(v,edir);

    // bring RA into conventional range
    edir.set_ra(dranrm(edir.get_ra()));
}

}
