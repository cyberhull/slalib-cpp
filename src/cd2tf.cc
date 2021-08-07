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
 * Converts an interval in days into hours, minutes, seconds (single precision).
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param ndp Number of decimal places of seconds; negative value interpreted as zero; for `days` up to 1.0, the
 *   available floating-point precision corresponds roughly to `ndp`=3.
 * @param days Interval in days; the absolute value of `days` may exceed 1.0; in cases where it does not, it is up to
 *   the caller to test for and handle the case where `days` is very nearly 1.0 and rounds up to 24 hours, by testing
 *   for result.get_hours()=24 and setting hours, minutes, seconds, and fraction to zero.
 * @param result Conversion result: hours, minutes, seconds, fraction, and sign
 */
void cd2tf(int ndp, float days, ConversionResult& result) {
    dd2tf(ndp, double (days), result);
}

}
