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
 * Transforms from J2000.0 equatorial coordinates to ecliptic coordinates (double precision).
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param dir J2000.0 mean RA,Dec (radians).
 * @param date TDB (barycentric dynamical time; loosely ET) as Modified Julian Date (JD-2400000.5).
 * @param edir Return value: ecliptic longitude and latitude; mean of date, IAU 1980 theory (radians).
 */
void eqecl(const Spherical<double>& dir, double date, Spherical<double>& edir) {
    Matrix<double> mat;
    Vector<double> v1, v2;

    // spherical to Cartesian
    dcs2c(dir, v1);

    // mean J2000 to mean of date
    prec(2000.0, epj(date), mat);
    dmxv(mat, v1, v2);

    // equatorial to ecliptic
    ecmat(date, mat);
    dmxv(mat, v2, v1);

    // Cartesian to spherical
    dcc2s(v1, edir);

    // express in conventional ranges
    edir.set_longitude(dranrm(edir.get_longitude()));
    edir.set_latitude(drange(edir.get_latitude()));
}

}
