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
 * Converts from universal time to sidereal time (double precision).
 *
 * The IAU 1982 expression (see page S15 of 1984 Astronomical Almanac) is used, but rearranged to reduce rounding
 * errors. This expression is always described as giving the GMST at 0 hours UT. In fact, it gives the difference
 * between the GMST and the UT, which happens to equal the GMST (modulo 24 hours) at 0 hours UT each day. In this
 * function, the entire UT is used directly as the argument for the standard formula, and the fractional part of the
 * UT is added separately. Note that the factor 1.0027379... does not appear in the IAU 1982 expression explicitly
 * but in the form of the coefficient 8640184.812866, which is 86400x36525x0.0027379...
 *
 * See also the function sla::gmsta(), which delivers better numerical precision by accepting the UT date and time as
 * separate arguments.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param ut1 Universal time (strictly UT1) expressed as modified Julian Date (JD-2400000.5).
 * @return Greenwich mean sidereal time (radians).
 */
double gmst(double ut1) {
    constexpr double D2PI = 6.283185307179586476925286766559;
    constexpr double S2R = 7.272205216643039903848711535369e-5;

    // Julian centuries from fundamental epoch J2000 to this UT
    const double jc = (ut1 - 51544.5) / 36525.0;

    // GMST at this UT
    return dranrm(std::fmod(ut1, 1.0) * D2PI +
        (24110.54841 +
        (8640184.812866 +
        (0.093104 - 6.2e-6 * jc) * jc) * jc) * S2R);
}

}
