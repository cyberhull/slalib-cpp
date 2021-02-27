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
 * Original FORTRAN code by P.T. Wallace.
 *
 * Reference: Explanatory Supplement to the Astronomical Almanac, ed P.K.Seidelmann,
 * University Science Books (1992), p604-606.
 *
 * @param year Year in Gregorian calendar, must be -4711 (i.e. 4712BC) or later.
 * @param month Month in Gregorian calendar; must be within [1..12] range.
 * @param day Day in Gregorian calendar; must be in the range [1..<days-in-given-month>]; if the day is out of range,
 *   `jdate` and `jday` are still computed and returned (but returned status is G2J_BAD_DAY).
 * @param jyear Output: year (re-aligned Julian calendar).
 * @param jday Output: day in year (1 = January 1st).
 * @return Conversion status as member of the G2JStatus enum; G2J_OK means conversion was successful.
 */
G2JStatus clyd(int year, int month, int day, int& jyear, int& jday) {
    // validate year
    if (year >= -4711) {
        // validate month
        if (month >= 1 && month <= 12) {
            // perform the conversion
            int i = (14 - month) / 12;
            int k = year - i;
            int j = (1461 * (k + 4800)) / 4 + (367 * (month - 2 + 12 * i)) / 12 -
                    (3 * ((k + 4900) / 100)) / 4 + day - 30660;
            k = (j - 1) / 1461;
            const int L = j - 1461 * k;
            const int N = (L - 1) / 365 - L / 1461;
            j = ((80 * (L - 365 * N + 30)) / 2447) / 11;
            i = N + j;
            jday = 59 + L - 365 * i + ((4 - N) / 4) * (1 - j);
            jyear = 4 * k + i - 4716;

            // validate day
            int days_in_given_month;
            if (month == 2) { // February?
                // allow for leap year
                days_in_given_month = (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0))? 29: 28;
            } else {
                // months' lengths in days
                static int days_in_month[12] = {
                    31, 0 /* February, handled above */, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
                };
                days_in_given_month = days_in_month[month - 1];
            }
            return (day < 1 || day > days_in_given_month)? G2J_BAD_DAY: G2J_OK;
        } else {
            return G2J_BAD_MONTH;
        }
    } else {
        return G2J_BAD_YEAR;
    }
}

}
