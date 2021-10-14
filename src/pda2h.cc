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
#include "f77_utils.h"
#include <cmath>

namespace sla {

/**
 * Calculates hour angle corresponding to a given azimuth (double precision).
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param lat Latitude (radians).
 * @param dec Declination (radians).
 * @param azimuth Azimuth (radians).
 * @param ha1 Return value: hour angle; first solution, if any (radians).
 * @param ha1_valid Return value: flag, `true` if first solution is valid (boolean).
 * @param ha2 Return value: hour angle; second solution, if any (radians).
 * @param ha2_valid Return value: flag, `true` if second solution is valid (boolean).
 */
void pda2h(double lat, double dec, double azimuth, double& ha1, bool& ha1_valid, double& ha2, bool& ha2_valid) {
    constexpr double PI = 3.141592653589793238462643;
    constexpr double DEG_90 = PI / 2.0;
    constexpr double TINY = 1.0e-12;

    // adjust latitude, azimuth, declination to avoid critical values
    double pn = drange(lat);
    if (std::fabs(std::fabs(pn) - DEG_90) < TINY) {
        pn -= f_sign(TINY, pn);
    } else if (std::fabs(pn) <TINY) {
        pn = TINY;
    }
    double an = drange(azimuth);
    if (std::fabs(std::fabs(an) - PI) < TINY) {
        an -= f_sign(TINY, an);
    } else if (std::fabs(an) < TINY) {
        an = TINY;
    }
    double dn = drange(dec);
    if (std::fabs(std::fabs(dn) - std::fabs(lat)) < TINY ||
        std::fabs(std::fabs(dn) - DEG_90) < TINY) {
        dn -= f_sign(TINY, dn);
    } else if (std::fabs(dn) < TINY) {
        dn = TINY;
    }

    // useful functions
    const double sin_an = std::sin(an);
    const double cos_an = std::cos(an);
    const double sinan_cospn = sin_an * std::sin(pn);

    // quotient giving sin(h+t)
    const double qt = std::sin(dn) * sin_an * std::cos(pn);
    const double qb = std::cos(dn) * std::sqrt(cos_an * cos_an + sinan_cospn * sinan_cospn);

    // any solutions?
    if (std::fabs(qt) <= qb) {

        // yes: find h+t and t
        const double hpt = std::asin(qt / qb);
        const double t = std::atan2(sinan_cospn, -cos_an);

        // the two solutions
        ha1 = drange(hpt - t);
        ha2 = drange(-hpt - (t + PI));

        // reject unless h and azimuth different signs
        ha1_valid = ha1 * an <= 0.0;
        ha2_valid = ha2 * an <= 0.0;
    } else {
        ha1_valid = ha2_valid = false;
    }
}

}
