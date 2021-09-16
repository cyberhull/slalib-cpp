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
 * Computes the E-terms (elliptic component of annual aberration) vector (double precision).
 *
 * Note the use of the J2000 aberration constant (20.49552 arcsec). This is a reflection of the fact that the E-terms
 * embodied in existing star catalogues were computed from a variety of aberration constants. Rather than adopting
 * one of the old constants the latest value is used here.
 *
 * References:
 *   Smith, C.A. et al., 1989.  Astr.J. 97, 265.
 *   Yallop, B.D. et al., 1989.  Astr.J. 97, 274.
 *
 * @param be Besselian epoch.
 * @param et Return value: E-terms as (dx, dy, dz).
 */
void etrms(double be, Vector<double> et) {
    // arcseconds to radians
    constexpr double ARCSECS_2_RADIANS = 0.484813681109535994e-5;

    // Julian centuries since B1950
    const double jc = (be - 1950.0) * 1.00002135903e-2;

    // eccentricity
    const double e = 0.01673011 - (0.00004193 + 0.000000126 * jc) * jc;

    // mean obliquity
    const double e0 = (84404.836 - (46.8495 + (0.00319 + 0.00181 * jc) * jc) * jc) * ARCSECS_2_RADIANS;

    // mean longitude of perihelion
    const double pl = (1015489.951 + (6190.67 + (1.65 + 0.012 * jc) * jc) * jc) * ARCSECS_2_RADIANS;

    // e-terms
    const double ek = e * 20.49552 * ARCSECS_2_RADIANS;
    const double cp = std::cos(pl);
    et[0] = ek * std::sin(pl);
    et[1] = -ek * cp * std::cos(e0);
    et[2] = -ek * cp * std::sin(e0);
}

}
