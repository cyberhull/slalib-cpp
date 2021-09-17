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
#include <cmath>

namespace sla {

/**
 * Converts from Universal Time to Greenwich mean sidereal time, with rounding errors minimized (double precision).
 *
 * There is no restriction on how the UT is apportioned between the `date` and `fdate` arguments. Either of the two
 * arguments could, for example, be zero and the entire date+time supplied in the other. However, the function is
 * designed to deliver maximum accuracy when the `date` argument is a whole number and the `fdate` lies in the range
 * 0 to 1 (or vice versa).
 *
 * The algorithm is based on the IAU 1982 expression (see page S15 of the 1984 Astronomical Almanac). This is always
 * described as giving the GMST at 0 hours UT1. In fact, it gives the difference between the GMST and the UT, the
 * steady 4-minutes-per-day drawing-ahead of ST with respect to UT. When whole days are ignored, the expression
 * happens to equal the GMST at 0 hours UT1 each day. Note that the factor 1.0027379... does not appear explicitly
 * but in the form of the coefficient 8640184.812866, which is 86400x36525x0.0027379...
 *
 * In this function, the entire UT1 (the sum of the two arguments `date` and `fdate`) is used directly as the argument
 * for the standard formula. The UT1 is then added, but omitting whole days to conserve accuracy.
 *
 * See also the function sla::gmst(), which accepts the UT as a single argument. Compared with sla::gmst(), the extra
 * numerical precision delivered by the present function is unlikely to be important in an absolute sense, but may be
 * useful when critically comparing algorithms and in applications where two sidereal times close together are
 * differenced.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param date UT1 date (MJD: integer part of JD-2400000.5).
 * @param fdate UT1 time (fraction of a day).
 * @return Greenwich mean sidereal time (radians, in the range 0..2pi).
 */
double gmsta(double date, double fdate) {
    // seconds of date to radians
    constexpr double SECONDS_2_RADIANS = 7.272205216643039903848712e-5;

    // Julian centuries since J2000
    double d1, d2;
    if (date < fdate) {
        d1 = date;
        d2 = fdate;
    } else {
        d1 = fdate;
        d2 = date;
    }
    const double jc = (d1 + (d2 - 51544.5)) / 36525.0;

    // GMST at this UT1
    return dranrm(SECONDS_2_RADIANS * (24110.54841 + (8640184.812866 +
        (0.093104 - 6.2e-6 * jc) * jc) * jc + 86400.0 * (std::fmod(d1, 1.0) + std::fmod(d2, 1.0))));
}

}
