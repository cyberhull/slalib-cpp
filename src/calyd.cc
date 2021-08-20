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
 * Converts Gregorian calendar to year and day in year (in a Julian calendar aligned to the 20th/21st century
 * Gregorian calendar).
 *
 * This routine exists to support the low-precision procedures sla::earth(), sla::moon(), and sla::ecor(). Between
 * 1900 March 1 and 2100 February 28 it returns answers which are consistent with the ordinary Gregorian calendar.
 * Outside this range there will be a discrepancy which increases by one day for every non-leap century year.
 *
 * The essence of the algorithm is first to express the Gregorian date as a Julian Day Number and then to convert
 * this back to a Julian calendar date, with day-in-year instead of month and day. See 12.92-1 and 12.95-1 in the
 * reference.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * Reference: Explanatory Supplement to the Astronomical Almanac, ed P.K.Seidelmann,
 * University Science Books (1992), p604-606.
 *
 * @param year Year in Gregorian calendar, must be -4711 (i.e. 4712BC) or later. some ranges get special treatment: if
 *   the year is in [0..49] then it is interpreted as 2000..2049; and if the year is in [50..99], it is interpreted
 *   as 1950..1999.
 * @param month Month in Gregorian calendar; must be within [1..12] range.
 * @param day Day in Gregorian calendar; must be in the range [1..<days-in-given-month>]; if the day is out of range,
 *   `jdate` and `jday` are still computed and returned (but returned status is G2J_BAD_DAY).
 * @param jyear Output: year (re-aligned Julian calendar).
 * @param jday Output: day in year (1 = January 1st).
 * @return Conversion status as member of the G2JStatus enum; G2J_OK means conversion was successful.
 */
G2JStatus calyd(int year, int month, int day, int& jyear, int& jday) {
    // do the conversion and return Julian year and day in the year in `jyear` and `jday`, respectively
    return clyd(process_year_defaults(year), month, day, jyear, jday);
}

}
