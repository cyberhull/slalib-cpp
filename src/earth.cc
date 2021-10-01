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
 * Calculates approximate heliocentric position and velocity of the Earth (single precision).
 *
 * The date and time is TDB (Barycentric Dynamical Time; loosely ET) in a Julian calendar, which has been aligned to
 * the ordinary Gregorian calendar for the interval 1900 March 1 to 2100 February 28. The year and day can be
 * obtained by calling sla::calyd() or sla::clyd() functions.
 *
 * Max/RMS errors 1950-2050:
 *   13/5 E-5 AU = 19200/7600 km in position,
 *   47/26 E-10 AU/s = 0.0070/0.0039 km/s in speed.
 *
 * More accurate results are obtainable with the functions sla::evp() and sla::epv().
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param year Year.
 * @param day Day in year (1 = January 1-st).
 * @param fraction Fraction of day.
 * @param pv Return value: Earth position and velocity vector; represents mean equator and equinox of date; position
 *   part is in AU; velocity part is in AU/sec.
 */
void earth(int year, int day, float fraction, VectorPV<float>& pv) {
    constexpr float TWO_PI = 6.28318530718f;

    // mean orbital speed of Earth, AU/s
    constexpr float orbital_speed = 1.9913e-7;

    // mean Earth:EMB (Earth-Moon barycenter) distance and speed, AU and AU/s
    constexpr float emb_dist = 3.12e-5;
    constexpr float emb_speed = 8.31e-11;

    // whole years & fraction of year, and years since J1900.0
    auto year_since_1900 = (const float) (year - 1900);
    const int y4 = ((year % 4) + 4) % 4;
    const float yf = ((float) (4 * (day - 1 / (y4 + 1)) - y4 - 2) + 4.0f * fraction) / 1461.0f;
    const float t = year_since_1900 + yf;

    // geometric mean longitude of Sun (cf 4.881627938+6.283319509911*t MOD 2PI)
    const float elm = std::fmod(4.881628f + TWO_PI * yf + 0.00013420f * t, TWO_PI);

    // mean longitude of perihelion
    const float gamma = 4.908230f + 3.0005e-4f * t;

    // mean anomaly
    const float em = elm - gamma;

    // mean obliquity
    const float eps0 = 0.40931975f - 2.27e-6f * t;

    // eccentricity
    const float e = 0.016751f - 4.2e-7f * t;
    const float e_squared = e * e;

    // true anomaly
    const float v = em + 2.0f * e * std::sin(em) + 1.25f * e_squared * std::sin(2.0f * em);

    // true ecliptic longitude
    const float elt = v + gamma;

    // true distance
    const float r = (1.0f - e_squared) / (1.0f + e * std::cos(v));

    // Moon's mean longitude
    const float elmm = std::fmod(4.72f + 83.9971f * t, TWO_PI);

    // useful functions
    const float cos_elt = std::cos(elt);
    const float sin_eps0 = std::sin(eps0);
    const float cos_eps0 = std::cos(eps0);
    const float w1 = -r * std::sin(elt);
    const float w2 = -orbital_speed * (cos_elt + e * std::cos(gamma));
    const float sin_elmm = std::sin(elmm);
    const float cos_elmm = std::cos(elmm);

    // Earth position and velocity
    pv.set_x(-r * cos_elt - emb_dist * cos_elmm);
    pv.set_y((w1 - emb_dist * sin_elmm) * cos_eps0);
    pv.set_z(w1 * sin_eps0);
    pv.set_dx(orbital_speed * (std::sin(elt) + e * std::sin(gamma)) + emb_speed * sin_elmm);
    pv.set_dy((w2 - emb_speed * cos_elmm) * cos_eps0);
    pv.set_dz(w2 * sin_eps0);
}

}
