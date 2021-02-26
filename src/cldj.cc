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
 * Converts Gregorian Calendar date to Modified Julian Date.
 *
 * Original FORTRAN code by P.T. Wallace. The algorithm is adapted from Hatcher 1984 (QJRAS 25, 53-55).
 *
 * @param year Gregorian calendar year, must be -4699 (i.e. 4700BC) or later.
 * @param month Month, must be in range [1..12].
 * @param day Day of the month, must be in range [1..<days-in-given-month>]; if the day is out of range, Modified
 *   Julian Date is still computed and returned (but returned status is G2J_BAD_DAY).
 * @param mjd Output: modified Julian Date (JD-2400000.5) for 0 hrs.
 * @return Conversion status as member of the G2JStatus enum; G2J_OK means conversion was successful.
 */
G2JStatus cldj(int year, int month, int day, double& mjd) {
    // validate year
    if (year >= -4699 ) {
        // validate month
        if (month >= 1 && month <= 12) {
            // calculate and return Modified Julian Date
            mjd = double((1461 * (year - (12 - month) / 10 + 4712)) / 4
                + (306 * ((month + 9) % 12) + 5) / 10
                - (3 * ((year - (12 - month) / 10 + 4900) / 100)) / 4
                + day - 2399904);
            // validate day
            int days_in_given_month;
            if (month == 2) { // February?
                // allow for leap year
                days_in_given_month = (year % 4 == 0)? 29: 28;
                if (year % 100 == 0 && year % 400 != 0) days_in_given_month = 28;
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
