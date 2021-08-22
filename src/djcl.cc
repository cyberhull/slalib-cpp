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
 * Converts Modified Julian Date to Gregorian year, month, day, and fraction of a day.
 *
 * The algorithm is adapted from Hatcher 1984 (QJRAS 25, 53-55).
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param mjd Modified Julian Date (JD-2400000.5).
 * @param date Return value: `Date` structure holding year, month [1..12], day [1..365], and fraction of the
 *   day [0..1.0f]
 * @return `false` if there was no error, `true` if `mjd` was out of range (before 4701BC March 1).
 */
bool djcl(double mjd, Date& date) {
    // check if date is acceptable.
    if (mjd <= -2395520.0 || mjd >= 1.0e9) {
        return true;
    } else {
        // separate day and fraction
        double f = std::fmod(mjd, 1.0);
        if (f < 0.0) {
            f += 1.0;
        }
        const double d = f_anint(mjd - f);

        // express day in Gregorian calendar.
        const int jd = f_nint(d) + 2400001;

        const int n4 = 4 * (jd + ((6 * ((4 * jd - 17918) / 146097)) / 4 + 1) / 2 - 37);
        const int nd10 = 10 * (((n4 - 237) % 1461) / 4) + 5;

        date.d_year = n4 / 1461 - 4712;
        date.d_month = ((nd10 / 306 + 2) % 12) + 1;
        date.d_day = (nd10 % 306) / 10 + 1;
        date.d_fraction = (float) f;
        return false;
    }
}

}
