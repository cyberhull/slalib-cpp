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
 * Returns increment to be applied to Coordinated Universal Time UTC to give International Atomic Time TAI (double
 * precision).
 *
 *   ========================================================================
 *   |                                                                      |
 *   |                              IMPORTANT!                              |
 *   |                                                                      |
 *   |  This function must be updated each time a leap second is announced  |
 *   |                                                                      |
 *   |                   Latest leap second: 2016 December 31               |
 *   |                                                                      |
 *   ========================================================================
 *
 * In original FORTRAN code, for epochs from 1961 January 1 onwards, the expressions from the file
 *   ftp://maia.usno.navy.mil/ser7/tai-utc.dat
 * were used; that link seem to be inaccessible as of now (August 15, 2021). The 5ms time step at 1961 January 1 is
 * taken from 2.58.1 (p87) of the 1992 Explanatory Supplement.
 *
 * C++ code updates may take data from the
 *   https://www.ietf.org/timezones/data/leap-seconds.list
 * file. The timestamps for each leap second introduction are at the bottom pf the file, and are in the NTP timestamp
 * format (seconds since the NTP epoch, which is 1 January 1900, 00:00:00). The `if` statements in the functions
 * compare `uts` argument to the Modified Julian Day number, which can be calculated from the NTP timestamp as follows:
 *
 *   <modified-julian-day-number> = <ntp-timestamp> / 86400 + 15020
 *
 * For instance, 3692217600 (NTP timestamp for January 1st 2017) yields MJD number 3692217600/86400+15020==57754.0,
 * which is the constant used in the (currently!) first `if` statement.
 *
 * No leap second will be introduced at the end of December 2021. See
 *   https://datacenter.iers.org/data/latestVersion/16_BULLETIN_C16.txt
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param utc UTC date as a modified JD (JD-2400000.5); the UTC is specified to be a date rather than a time to
 *   indicate that care needs to be taken not to specify an instant which lies within a leap second; though in most
 *   cases UTC can include the fractional part, correct behaviour on the day of a leap second can only be guaranteed
 *   up to the end of the second 23:59:59; UTC began at 1960 January 1.0 (JD 2436934.5) and it is improper to call
 *   this function with an earlier epoch; however, if this is attempted, the TAI-UTC expression for the year 1960 is
 *   used.
 * @return TAI minus UTC (seconds).
 */
double dat(double utc) {

    /***********************************************/
    /*  Add new code here on each occasion that a  */
    /*  leap second is announced, and update the   */
    /*  preamble comments appropriately.           */
    /***********************************************/

    if (utc >= 57754.0) { // 2017 January 1
        return 37.0;
    } else if (utc >= 57204.0) { // 2015 July 1
        return 36.0;
    } else if (utc >= 56109.0) { // 2012 July 1
        return 35.0;
    } else if (utc >= 54832.0) { // 2009 January 1
        return 34.0;
    } else if (utc >= 53736.0) { // 2006 January 1
        return 33.0;
    } else if (utc >= 51179.0) { // 1999 January 1
        return 32.0;
    } else if (utc >= 50630.0) { // 1997 July 1
        return 31.0;
    } else if (utc >= 50083.0) { // 1996 January 1
        return 30.0;
    } else if (utc >= 49534.0) { // 1994 July 1
        return 29.0;
    } else if (utc >= 49169.0) { // 1993 July 1
        return 28.0;
    } else if (utc >= 48804.0) { // 1992 July 1
        return 27.0;
    } else if (utc >= 48257.0) { // 1991 January 1
        return 26.0;
    } else if (utc >= 47892.0) { // 1990 January 1
        return 25.0;
    } else if (utc >= 47161.0) { // 1988 January 1
        return 24.0;
    } else if (utc >= 46247.0) { // 1985 July 1
        return 23.0;
    } else if (utc >= 45516.0) { // 1983 July 1
        return 22.0;
    } else if (utc >= 45151.0) { // 1982 July 1
        return 21.0;
    } else if (utc >= 44786.0) { // 1981 July 1
        return 20.0;
    } else if (utc >= 44239.0) { // 1980 January 1
        return 19.0;
    } else if (utc >= 43874.0) { // 1979 January 1
        return 18.0;
    } else if (utc >= 43509.0) { // 1978 January 1
        return 17.0;
    } else if (utc >= 43144.0) { // 1977 January 1
        return 16.0;
    } else if (utc >= 42778.0) { // 1976 January 1
        return 15.0;
    } else if (utc >= 42413.0) { // 1975 January 1
        return 14.0;
    } else if (utc >= 42048.0) { // 1974 January 1
        return 13.0;
    } else if (utc >= 41683.0) { // 1973 January 1
        return 12.0;
    } else if (utc >= 41499.0) { // 1972 July 1
        return 11.0;
    } else if (utc >= 41317.0) { // 1972 January 1
        return 10.0;
    } else if (utc >= 39887.0) { // 1968 February 1
        return 4.21317 + (utc - 39126.0) * 0.002592;
    } else if (utc >= 39126.0) { // 1966 January 1
        return 4.31317 + (utc - 39126.0) * 0.002592;
    } else if (utc >= 39004.0) { // 1965 September 1
        return 3.84013 + (utc - 38761.0) * 0.001296;
    } else if (utc >= 38942.0) { // 1965 July 1
        return 3.74013 + (utc - 38761.0) * 0.001296;
    } else if (utc >= 38820.0) { // 1965 March 1
        return 3.64013 + (utc - 38761.0) * 0.001296;
    } else if (utc >= 38761.0) { // 1965 January 1
        return 3.54013 + (utc - 38761.0) * 0.001296;
    } else if (utc >= 38639.0) { // 1964 September 1
        return 3.44013 + (utc - 38761.0) * 0.001296;
    } else if (utc >= 38486.0) { // 1964 April 1
        return 3.34013 + (utc - 38761.0) * 0.001296;
    } else if (utc >= 38395.0) { // 1964 January 1
        return 3.24013 + (utc - 38761.0) * 0.001296;
    } else if (utc >= 38334.0) { // 1963 November 1
        return 1.945858 + (utc - 37665.0) * 0.0011232;
    } else if (utc >= 37665.0) { // 1962 January 1
        return 1.845858 + (utc - 37665.0) * 0.0011232;
    } else if (utc >= 37512.0) { // 1961 August 1
        return 1.372818 + (utc - 37300.0) * 0.001296;
    } else if (utc >= 37300.0) { // 1961 January 1
        return 1.422818 + (utc - 37300.0) * 0.001296;
    } else { // Before that
        return 1.417818 + (utc - 37300.0) * 0.001296;
    }
}

}
