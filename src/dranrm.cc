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
 * Normalize angle into range 0-2*pi (double precision).
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param angle The angle in radians.
 * @return The `angle` expressed in the range 0-2*pi .
 */
double dranrm(double angle) {
    constexpr double A2PI = 6.283185307179586476925287;
    double result = std::fmod(angle, A2PI);
    if (result < 0.0) {
        result += A2PI;
    }
    return result;
}

}
