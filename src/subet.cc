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
 * Removee the E-terms (elliptic component of annual aberration) from a pre IAU 1976 catalogue RA,Dec to give a mean
 * place (double precision).
 *
 * Most star positions from pre-1984 optical catalogues (or derived from astrometry using such stars) embody the
 * E-terms. This routine converts such a position to a formal mean place (allowing, for example, comparison with a
 * pulsar timing position).
 *
 * Reference:
 *   Explanatory Supplement to the Astronomical Ephemeris, section 2D, page 48.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param edir RA,Dec with E-terms included (radians).
 * @param be Besselian epoch of mean equator and equinox.
 * @param dir Return value: RA,Dec without E-terms (radians).
 */
void subet(const Spherical<double>& edir, double be, Spherical<double>& dir) {
    // get E-terms
    Vector<double> et;
    etrms(be, et);

    // spherical to Cartesian
    Vector<double> v;
    dcs2c(edir, v);

    // include the E-terms
    const double f = 1.0 + dvdv(v, et);
    for (int i = 0; i < 3; i++) {
        v[i] = v[i] * f - et[i];
    }
    // Cartesian to spherical
    dcc2s(v, dir);

    // bring RA into conventional range
    dir.set_ra(dranrm(dir.get_ra()));
}

}
