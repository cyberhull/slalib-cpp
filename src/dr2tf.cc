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
 * Convert an angle in radians to hours, minutes, seconds (double precision).
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param ndp Number of decimal places of seconds; negative value interpreted as zero; for `angle` up to 2*Pi, the
 *   available floating-point precision corresponds roughly to `ndp`==12.
 * @param angle Angle in radians; the absolute value of `angle` may exceed 2*Pi; in cases where it does not, it is up
 *   to the caller to test for and handle the case where `angle` is very nearly 2*Pi and rounds up to 24 hours, by
 *   testing for result.get_hours()==24 and setting hours, minutes, seconds, and fraction to zero.
 * @param result Conversion result: hours, minutes, seconds, fraction, and sign.
 */
void dr2tf(int ndp, double angle, ConversionResult& result) {
    // turns to radians
    constexpr double TURNS2RADIANS = 6.283185307179586476925287;
    // scale, then use days to hours, minutes, seconds function
    dd2tf(ndp, angle / TURNS2RADIANS, result);
}

}
