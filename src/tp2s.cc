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
 * Transforms tangent plane coordinates into spherical (single precision).
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param xi Tangent plane rectangular coordinate.
 * @param eta Tangent plane rectangular coordinate.
 * @param tangent Spherical coordinates of tangent point.
 * @param point Return value: spherical coordinates (0-2pi,+/-pi/2).
 */
void tp2s(float xi, float eta, const Spherical<float>& tangent, Spherical<float>& point) {
    const float sin_tdec = std::sin(tangent.get_dec());
    const float cos_tdec = std::cos(tangent.get_dec());
    const float denom = cos_tdec - eta * sin_tdec;

    point.set_ra(ranorm(std::atan2(xi, denom) + tangent.get_ra()));
    point.set_dec(std::atan2(sin_tdec + eta * cos_tdec, std::sqrt(xi * xi + denom * denom)));
}

}
