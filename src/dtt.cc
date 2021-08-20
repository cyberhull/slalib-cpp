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
 * Returns increment to be applied to Coordinated Universal Time UTC to give Terrestrial Time TT (formerly Ephemeris
 * Time ET) (double precision).
 *
 * See also the sla::dt() function, which roughly estimates ET-UT for historical epochs.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param utc UTC date as a modified JD (JD-2400000.5); The UTC is specified to be a date rather than a time to
 *   indicate that care needs to be taken not to specify an instant which lies within a leap second. Though in most
 *   cases UTC can include the fractional part, correct behaviour on the day of a leap second can only be guaranteed
 *   up to the end of the second 23:59:59.
 * @return TT minus UTC (seconds); pre 1972 January 1 a fixed value of 10 + ET-TAI is returned.
 */
double dtt(double utc) {
    return 32.184 + dat(utc);
}

}
