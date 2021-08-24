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
 * Calculates velocity component in a given direction due to Earth rotation (single precision).
 *
 * The simple algorithm used assumes a spherical Earth, of a radius chosen to give results accurate to about 0.0005
 * km/s for observing stations at typical latitudes and heights. For applications requiring greater precision, use the
 * function sla::pvobs().
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param phi Latitude of observing station (geodetic; radians).
 * @param pos Apparent RA,Dec (radians).
 * @param stime Local apparent sidereal time (radians).
 * @return Component of Earth rotation in direction RA,DA (km/s); the result is positive when the observatory is
 *   receding from the given point on the sky.
 */
float rverot(float phi, const Spherical<float>& pos, float stime) {
    // nominal mean sidereal speed of Earth equator in km/s (the actual value is about 0.4651).
    constexpr float EARTH_SPEED = 0.4655f;
    return EARTH_SPEED * std::cos(phi) * std::sin(stime - pos.get_ra()) * std::cos(pos.get_dec());
}

}
