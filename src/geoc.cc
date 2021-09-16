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
 * Converts geodetic position to geocentric (double precision).
 *
 * Uses IAU 1976 constants.
 *
 * After getting return values, geocentric latitude can be calculated as `atan2(equator_dist, axis_dist)`.
 *
 * Reference:
 *   Green,R.M., Spherical Astronomy, CUP 1985, p98.
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param latitude Geodetic latitude (radians).
 * @param height Geodetic height above reference spheroid (meters).
 * @param axis_dist Return value: distance from Earth axis (AU).
 * @param equator_dist Return value: distance from plane of Earth equator (AU).
 */
void geoc(double latitude, double height, double& axis_dist, double& equator_dist) {
    // Earth equatorial radius (meters)
    constexpr double EARTH_RADIUS = 6378140.0;

    // reference spheroid flattening factor and useful function
    constexpr double FF = 1.0 / 298.257;
    constexpr double B2 = 1.0 - FF;
    constexpr double B = B2 * B2;

    // astronomical unit in meters
    constexpr double AU = 1.49597870e11;

    // geodetic to geocentric conversion
    const double sin_lat = std::sin(latitude);
    const double cos_lat = std::cos(latitude);
    const double c = 1.0 / std::sqrt(cos_lat * cos_lat + B * sin_lat * sin_lat);
    const double s = B * c;

    // return result
    axis_dist = (EARTH_RADIUS * c + height) * cos_lat / AU;
    equator_dist = (EARTH_RADIUS * s + height) * sin_lat / AU;
}

}
