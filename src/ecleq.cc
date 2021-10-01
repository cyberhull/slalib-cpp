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
 * Transforms ecliptic coordinates to J2000.0 equatorial coordinates (double precision).
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param ecliptic Ecliptic longitude and latitude (mean of date, IAU 1980 theory, radians).
 * @param date TDB (Barycentric Dynamical Time; loosely ET) as Modified Julian Date (JD - 2400000.5).
 * @param equatorial J2000.0 mean RA,Dec (radians).
 */
void ecleq(const Spherical<double>& ecliptic, double date, Spherical<double>& equatorial) {
    // spherical to Cartesian
    Vector<double> v1;
    dcs2c(ecliptic, v1);

    // ecliptic to equatorial
    Matrix<double> mat;
    ecmat(date, mat);
    Vector<double> v2;
    dimxv(mat, v1, v2);

    // mean of date to J2000
    prec(2000.0, epj(date), mat);
    dimxv(mat, v2, v1);

    // Cartesian to spherical
    dcc2s(v1, equatorial);

    // express in conventional ranges
    equatorial.set_ra(dranrm(equatorial.get_ra()));
    equatorial.set_dec(drange(equatorial.get_dec()));
}

}
