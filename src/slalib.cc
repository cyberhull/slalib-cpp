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

int process_year_defaults(int year) {
    if (year >= 0 && year <= 49) {
        return year + 2000;
    } else if (year >= 50 && year <= 99) {
        return year + 1900;
    }
    return year;
}

G2JStatus validate_gregorian_day(int year, int month, int day) {
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
}

}
