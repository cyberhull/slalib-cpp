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

namespace sla {

/**
 * Normalizes angle into range +/- pi  (single precision).
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param angle The angle in radians.
 * @return The `angle` expressed in the +/- pi.
 */
float range(float angle) {
    constexpr float PI = 3.141592653589793238462643f;
    constexpr float PI2 = 6.283185307179586476925287f;
    float result = std::fmod(angle, PI2);
    if (std::fabs(result) >= PI) {
        result -= f_sign(PI2, angle);
    }
    return result;
}

}
