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
 * Calculates and returns position and velocity of an observing station (double precision).
 *
 * Uses IAU 1976 constants.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param latitude Geodetic latitude (radians).
 * @param height Geodetic height above reference spheroid (meters).
 * @param lst Local apparent sidereal time (radians).
 * @param pv Return value: position/velocity vector (AU, AU/s, true equator and equinox of date).
 */
void pvobs(double latitude, double height, double lst, VectorPV<double>& pv) {
    // mean sidereal rate (at J2000) in radians per (UT1) second
    constexpr double SIDERIAL_RATE = 7.292115855306589e-5;

    // geodetic to geocentric conversion
    double axis_dist, equator_dist;
    geoc(latitude, height, axis_dist, equator_dist);

    // functions of lst
    const double sin_lst = std::sin(lst);
    const double cos_lst = std::cos(lst);

    // speed
    const double velocity = SIDERIAL_RATE * axis_dist;

    // position
    pv.set_x(axis_dist * cos_lst);
    pv.set_y(axis_dist * sin_lst);
    pv.set_z(equator_dist);

    // velocity
    pv.set_dx(-velocity * sin_lst);
    pv.set_dy(velocity * cos_lst);
    pv.set_dz(0.0);
}

}
