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
 * Forms the matrix of precession between two epochs (IAU 1976, FK5) (double precision).
 *
 *  Though the matrix method itself is rigorous, the precession angles are expressed through canonical polynomials
 *  which are valid only for a limited time span. There are also known errors in the IAU precession rate. The
 *  absolute accuracy of the present formulation is better than 0.1 arcsecond from 1960AD to 2040AD, better than
 *  1 arcsecond from 1640AD to 2360AD, and remains below 3 arcseconds for the whole of the period 500BC to 3000AD. The
 *  errors exceed 10 arcseconds outside the range 1200BC to 3900AD, exceed 100 arcseconds outside 4200BC to 5600AD and
 *  exceed 1000 arcseconds outside 6800BC to 8200AD. The SLALIB function sla::precl() implements a more elaborate
 *  model, which is suitable for problems spanning several thousand years.
 *
 * References:
 *    Lieske,J.H., 1979. Astron.Astrophys.,73,282.
 *     equations (6) & (7), p283.
 *    Kaplan,G.H., 1981. USNO circular no. 163, pA2.
 *
 * @param ep0 Beginning epoch: TDB (Barycentric Dynamical Time; loosely ET, Ephemeris Time) Julian epoch.
 * @param ep1 Ending epoch: TDB (Barycentric Dynamical Time; loosely ET, Ephemeris Time) Julian epoch.
 * @param mat Return value: precession matrix; the matrix is in the sense v(ep1) = mat * v(ep0).
 */
void prec(double ep0, double ep1, matrix<double> mat) {
    // arc seconds to radians
    constexpr double ARCSECS_2_RADIANS = 0.484813681109535994e-5;

    // interval between basic epoch J2000.0 and beginning epoch (JC)
    const double t0 = (ep0 - 2000.0) / 100.0;

    // interval over which precession required (JC)
    const double t = (ep1 - ep0) / 100.0;

    // Euler angles
    const double tas2R = t * ARCSECS_2_RADIANS;
    const double w = 2306.2181 + (1.39656 - 0.000139 * t0) * t0;

    const double zeta = (w + ((0.30188 - 0.000344 * t0) + 0.017998 * t) * t) * tas2R;
    const double z = (w + ((1.09468 + 0.000066 * t0) + 0.018203 * t) * t) * tas2R;
    const double theta = ((2004.3109 + (-0.85330 - 0.000217 * t0) * t0) + ((-0.42665 - 0.000217 * t0) -
        0.041833 * t) * t) * tas2R;

    // rotation matrix
    deuler("ZYZ", -zeta, theta, -z, mat);
}

}
