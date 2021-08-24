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
 * Forms the matrix of nutation for a given tdb - Shirai & Fukushima 2001 theory (double precision).
 *
 * Earth attitude predictions made by combining the present nutation matrix with IAU~1976 precession are accurate
 * to 1~milliarcsecond (with respect to the ICRS) for a few decades around 2000.
 *
 * The distinction between the required TDB (Barycentric Dynamical Time) and TT (Terrestrial Time) is always
 * negligible. Moreover, for all but the most critical applications UTC is adequate.
 *
 * Reference:
 *     Shirai, T. & Fukushima, T., Astron.J. 121, 3270-3283 (2001).
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param tdb TDB (Barycentric Dynamical Time; loosely ET, Ephemeris Time) as Modified Julian Date (JD-2400000.5).
 * @param mat Nutation matrix; the matrix is in the sense  v(true) = mat * v(mean), where v(true) is the star vector
 *   relative to the true equator and equinox of tdb and v(mean) is the star vector relative to the mean equator
 *   and equinox of tdb; the matrix represents forced nutation (but not free core nutation) plus corrections to the
 *   IAU~1976 precession model.
 */
void nut(double tdb, Matrix<double> mat) {
    double psi, eps, eps0;

    // nutation components and mean obliquity
    nutc(tdb,psi,eps,eps0);

    // rotation matrix
    deuler("XZX", eps0, -psi,-(eps0 + eps), mat);
}

}
