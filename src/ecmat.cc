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
 * Form the equatorial to ecliptic rotation matrix - IAU 1980 theory (double precision).
 *
 * Reference: Murray,C.A., Vectorial Astrometry, section 4.3.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param date (Loosely ET) as Modified Julian Date (JD-2400000.5).
 * @param mat Return value: rotation matrix; the matrix is in the sense  V(ecl) == RMAT * V(equ); the equator,
 *   equinox and ecliptic are mean of date.
 */
void ecmat(double date, Matrix<double> mat) {
    // conversion constant: arc seconds to radians
    constexpr double ARCSECONDS2RADIANS = 0.484813681109535994e-5;

    // interval between basic epoch J2000.0 and current epoch (JC)
    const double t = (date - 51544.5) / 36525.0;

    // mean obliquity
    const double phi = ARCSECONDS2RADIANS * (84381.448 + (-46.8150 + (-0.00059 + 0.001813 * t ) * t) * t);

    // calculate the matrix
    deuler("X", phi, 0.0, 0.0, mat);
}

}
