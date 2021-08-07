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
 * Convert hours, minutes, seconds to radians (single precision).
 *
 * The result is computed even if any of the range checks fail. The sign must be dealt with outside this routine.
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param hours Hours; must be in range [0..23].
 * @param minutes Minutes; must be in range [0..60].
 * @param seconds Seconds; must be in range [0..60).
 * @param radians Return value: angle in radians.
 * @return Conversion status: a T2D_xxx constant.
 */
T2DStatus ctf2r(int hours, int minutes, float seconds, float & radians) {
    T2DStatus result = ctf2d(hours, minutes, seconds, radians);
    constexpr float TURNS2RADIANS = 6.283185307179586476925287f;
    radians *= TURNS2RADIANS;
    return result;
}

}
