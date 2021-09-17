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
 * Calculates term for the equation of the equinoxes [IAU 1994] (double precision).
 *
 * References: IAU Resolution C7, Recommendation 3 (1994)
 *   Capitaine, N. & Gontier, A.-M., Astron. Astrophys., 275, 645-650 (1993).
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param date TDB (barycentric dynamical time; loosely ET) as Modified Julian Date (JD-2400000.5).
 * @return Term for the equation of the equinoxes (radians): Greenwich apparent ST = GMST + sla::eqeqx().
 */
double eqeqx(double date) {
    constexpr double TURNS_2_ARCSECS = 1296000.0;
    constexpr double ARCSECS_2_RADIANS = 0.484813681109535994e-5;

    // interval between basic epoch J2000.0 and current epoch (JC)
    const double t = (date - 51544.5) / 36525.0;

    // longitude of the mean ascending node of the lunar orbit on the
    // ecliptic, measured from the mean equinox of date
    const double om = ARCSECS_2_RADIANS *
        (450160.280 + (-5.0 * TURNS_2_ARCSECS - 482890.539 + (7.455 + 0.008 * t) * t) * t);

    // Nutation
    double psi, eps, eps0;
    nutc(date, psi, eps, eps0);

    // equation of the equinoxes
    return psi * std::cos(eps0) + ARCSECS_2_RADIANS * (0.00264 * std::sin(om) + 0.000063 * std::sin(om + om));
}

}
