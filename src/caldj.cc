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
 * Converts Gregorian Calendar date to Modified Julian Date, with support for century default. See `year` argument
 * description for details. Calls `sla::cldj()` to do the conversion. If the support for the default century is not
 * needed, or if one is converting a year before 100AD, `sla::cldj()` should be called directly.
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param year Gregorian calendar year, must be -4699 (i.e. 4700BC) or later; some ranges get special treatment: if
 *   the year is in [0..49] then it is interpreted as 2000..2049; and if the year is in [50..99], it is interpreted
 *   as 1950..1999.
 * @param month Month, must be in range [1..12].
 * @param day Day of the month, must be in range [1..<days-in-given-month>]; if the day is out of range, Modified
 *   Julian Date is still computed and returned (but returned status is G2J_BAD_DAY).
 * @param mjd Output: modified Julian Date (JD-2400000.5) for 0 hrs.
 * @return Conversion status as member of the G2JStatus enum; G2J_OK means conversion was successful.
 */
G2JStatus caldj(int year, int month, int day, double& mjd) {
    // do the conversion and return modified Julian Date in `mjd`
    return cldj(process_year_defaults(year), month, day, mjd);
}

}
