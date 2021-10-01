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
 * Calculates component of Earth orbit velocity and heliocentric light time in a given direction (single precision).
 *
 * The date and time is TDB (Barycentric Dynamical Time; loosely ET) in a Julian calendar which has been aligned to
 * the ordinary Gregorian calendar for the interval 1900 March 1 to 2100 February 28. The year and day can be
 * obtained by calling sla::calyd() or sla::clyd() functions.
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param dir Mean RA, Dec of date (radians).
 * @param year Year.
 * @param day Day in the year (1 = January 1-st).
 * @param fraction Fraction of day.
 * @param velocity Return value: component of Earth orbital velocity (km/sec); positive when the Earth is receding
 *   from the given point on the sky; accuracy is usually within 0.004 km/s of the correct value and is never in
 *   error by more than 0.007 km/s.
 * @param lt Return value: component of heliocentric light time (sec); positive when the Earth lies between the Sun
 *   and the given point on the sky; the error in light time correction is about 0.03s at worst, but is usually
 *   better than 0.01s; for applications requiring higher accuracy, see the sla::evp() and sla::epv() functions.
 */
void ecor(Spherical<float> dir, int year, int day, float fraction, float& velocity, float& lt) {
    // AU to km and light sec (1985 Almanac)
    constexpr float AU_2_KM = 1.4959787066e8f;
    constexpr float AU_2_LSEC = 499.0047837f;

    // Sun:Earth position & velocity vector
    VectorPV<float> pv;
    earth(year, day, fraction, pv);

    // star position vector
    Vector<float> vec;
    cs2c(dir, vec);

    // velocity component
    velocity = -AU_2_KM * vdv(pv.get_velocity(), vec);

    // light time component
    lt = AU_2_LSEC * vdv(pv.get_position(), vec);
}

}
