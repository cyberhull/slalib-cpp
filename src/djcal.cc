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
 * Converts Modified Julian Date to Gregorian Calendar, expressed in a form convenient for formatting messages (namely
 * rounded to a specified precision, and with the fields stored in a single structure).
 *
 * The algorithm is adapted from Hatcher 1984 (QJRAS 25, 53-55).
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param ndp Number of decimal places of days in fraction; should be 4 or less if internal overflows are to be avoided
 *   on machines which use 32-bit integers.
 * @param mjd Modified Julian Date (JD-2400000.5); any date after 4701BC March 1 is accepted.
 * @param date Return value: `Date` structure holding year, month, day, fraction (in the `d_ifraction` field) in
 *   Gregorian calendar.
 * @return Status: `true` means input date is out of range.
 */
bool djcal(int ndp, double mjd, Date& date) {
    // validate
    if (mjd <= -2395520.0 || mjd >= 1.0e9) {
        return true;
    } else {
        // denominator of fraction.
        int nfd = (int) std::round(std::pow(10.0, (double) std::max(ndp, 0)));
        const double fd = (double) nfd;

        // round date and express in units of fraction.
        const double df = f_anint(mjd * fd);

        // separate day and fraction.
        double f = std::fmod(df, fd);
        if (f < 0.0) {
            f += fd;
        }
        const double d = (df - f) / fd;

        // express day in Gregorian calendar.
        const int jd = f_nint(d) + 2400001;

        const int n4 = 4 * (jd + ((2 * ((4 * jd - 17918) / 146097) * 3) / 4 + 1) / 2 - 37);
        const int nd10 = 10 * (((n4 - 237) % 1461) / 4) + 5;

        date.d_year = n4 / 1461 - 4712;
        date.d_month =  ((nd10 / 306 + 2) % 12) + 1;
        date.d_day = (nd10 % 306) / 10 + 1;
        date.d_ifraction = f_nint(f);
        return false;
    }
}

}
