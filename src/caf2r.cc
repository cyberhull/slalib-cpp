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
 * Converts degrees, arc minutes, arc seconds to radians (single precision).
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param degrees Degrees in range [0..359]
 * @param minutes Arc minutes in range [0..59]
 * @param seconds Arc seconds in range [0..60.0)
 * @param radians Angle in radians; calculated and returned even if parameters are out of range.
 * @return Conversion status as member of the D2RStatus enum; D2R_OK means conversion was successful.
 */
D2RStatus caf2r(int degrees, int minutes, float seconds, float& radians) {
    // assume all arguments are within their ranges
    D2RStatus status = D2R_OK;

    // validate arguments
    if (seconds < 0.0f || seconds >= 60.0f) {
        status = D2R_BAD_ARCSECONDS;
    }
    if (minutes < 0 || minutes > 59) {
        status = D2R_BAD_ARCMINUTES;
    }
    if (degrees < 0 || degrees > 359) {
        status = D2R_BAD_DEGREES;
    }
    // calculate and return angle in radians
    constexpr float ARCSECS_2_RADIANS = 0.484813681109535994e-5f;
    radians = ((float(degrees) * 60.0f + float(minutes)) * 60.0f + seconds) * ARCSECS_2_RADIANS;
    return status;
}

}
