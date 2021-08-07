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
 * Convert hours, minutes, seconds to days (single precision).
 *
 * The result is computed even if any of the range checks fail. The sign must be dealt with outside this routine.
 *
 * Original FORTRAN code by Rutherford Appleton Laboratory / P.T. Wallace.
 *
 * @param hours Hours; must be in range [0..23].
 * @param minutes Minutes; must be in range [0..60].
 * @param seconds Seconds; must be in range [0..60).
 * @param days Return value: interval in days.
 * @return Conversion status: a T2D_xxx constant.
 */
T2DStatus ctf2d(int hours, int minutes, float seconds, float& days) {
    // validate arguments
    T2DStatus result = T2D_OK;
    if (hours < 0 || hours > 23) result = T2D_BAD_HOURS;
    if (minutes < 0 || minutes > 59) result = T2D_BAD_MINUTES;
    if (seconds < 0.0f || seconds >= 60.0f) result = T2D_BAD_SECONDS;

    // compute interval
    constexpr float DAYS2SECONDS = 86400.0f;
    days= (60.0f * (60.0f * (float) hours + (float) minutes) + seconds) / DAYS2SECONDS;
    return result;
}

}
