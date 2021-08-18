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
 * Forms the matrix of precession and nutation (SF2001) (double precision).
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param epoch Julian Epoch for mean coordinates; TDB (Barycentric Dynamical Time; loosely ET, Ephemeris Time); TT
 *   (Terrestrial Time) will also do, or even UTC.
 * @param date Modified Julian Date (JD-2400000.5) for true coordinates; TDB (Barycentric Dynamical Time; loosely ET,
 *   Ephemeris Time); TT (Terrestrial Time) will also do, or even UTC.
 * @param mat Return value: combined precession/nutation matrix; The matrix is in the sense v(true) = mat * v(mean).
 */
void prenut(double epoch, double date, matrix<double> mat) {
    matrix<double> pmat, nmat;

    // precession
    prec(epoch, epj(date), pmat);

    // nutation
    nut(date, nmat);

    // combine the matrices: PN = N x P
    dmxm(nmat, pmat, mat);
}

}
