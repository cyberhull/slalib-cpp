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
#include <cmath>

namespace sla {

/**
 * Converts an interval in days into hours, minutes, seconds (double precision).
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param ndp Number of decimal places of seconds; negative value interpreted as zero; for `days` up to 1.0, the
 *   available floating-point precision corresponds roughly to `ndp`=12.
 * @param days Interval in days; the absolute value of `days` may exceed 1.0; in cases where it does not, it is up to
 *   the caller to test for and handle the case where `days` is very nearly 1.0 and rounds up to 24 hours, by testing
 *   for result.get_hours()==24 and setting hours, minutes, seconds, and fraction to zero.
 * @param result Conversion result: hours, minutes, seconds, fraction, and sign
 */
void dd2tf(int ndp, double days, ConversionResult& result) {
    // field units in terms of least significant figure
    int nrs = 1;
    for (int i = 0; i < ndp; i++) {
        nrs *= 10;
    }
    auto rs = (double) nrs;
    double rm = rs * 60.0;
    double rh = rm * 60.0;

    // round interval and express in smallest units required
    constexpr double DAYS2SECONDS = 86400.0;
    double interval = f_anint(rs * DAYS2SECONDS * std::abs(days));

    // separate into fields
    double hours = f_aint(interval / rh);
    interval -= hours * rh;
    double minutes = f_aint(interval / rm);
    interval -= minutes * rm;
    double seconds = f_aint(interval / rs);
    double fraction = interval - seconds * rs;

    // return results
    result.set_hours(std::max(f_nint(hours), 0));
    result.set_minutes(std::max(std::min(f_nint(minutes), 59), 0));
    result.set_seconds(std::max(std::min(f_nint(seconds), 59), 0));
    result.set_fraction(std::max(f_nint(std::min(fraction, rs - 1.0)), 0));
    result.set_sign(days >= 0.0);
}

}
