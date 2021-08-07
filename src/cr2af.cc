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
 * Converts an angle in radians into degrees, arcminutes, arcseconds (single precision).
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param ndp Number of decimal places of arcseconds; negative value interpreted as zero; for `angle` up to 2*Pi, the
 *   available floating-point precision corresponds roughly to `ndp`==3.
 * @param angle Angle in radians; the absolute value of `angle` may exceed 2*Pi; in cases where it does not, it is up
 *   to the caller to test for and handle the case where `angle` is very nearly 2*Pi and rounds up to 360 degrees, by
 *   testing for result.degrees()==360 and setting degrees, arcminutes, arcseconds, and fraction to zero.
 * @param result Conversion result: degrees, arcminutes, arcseconds, fraction, and sign.
 */
void cr2af(int ndp, float angle, ConversionResult& result) {
    // hours to degrees * radians to turns
    constexpr auto factor = (float)(15.0 / 6.283185307179586476925287);
    cd2tf(ndp, angle * factor, result);
}

}
