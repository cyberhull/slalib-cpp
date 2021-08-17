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
 * Computes nutation: longitude & obliquity components and mean obliquity, using the IAU 1980 theory (double precision).
 *
 * Earth attitude predictions made by combining the present nutation model with IAU 1976 precession are accurate to
 * 0.35 arcseconds over the interval 1900-2100 (the accuracy is much better near the middle of the interval).
 *
 * The sla::nutc() function is the equivalent of the present function but uses the Shirai & Fukushima 2001 forced
 * nutation theory (SF2001). The newer theory is more accurate than IAU 1980, within 1 milliarcsecond (with respect to
 * the ICRF) for a few decades around 2000. The improvement is mainly because of the corrections to the IAU 1976
 * precession that the SF2001 theory includes.
 *
 * References:
 *   Final report of the IAU Working Group on Nutation, chairman P.K.Seidelmann, 1980.
 *   Kaplan,G.H., 1981, USNO circular no. 163, pA3-6.
 *
 * Original FORTRAN code by Rutherford Appleton Laboratory / P.T. Wallace.
 *
 * @param date TDB (Barycentric Dynamical Time; loosely ET, Ephemeris Time) as Modified Julian Date (JD-2400000.5).
 * @param psi Return value: nutation in longitude.
 * @param eps Return value: nutation in obliquity.
 * @param eps0 Return value: mean obliquity.
 */
void nutc80(double date, double& psi, double& eps, double& eps0) {
    // turns to arc seconds
    constexpr double TURNS_2_ARCSECS = 1296000.0;
    // arc seconds to radians
    constexpr double ARCSECS_2_RADIANS = 0.484813681109535994e-5;
    // units of 0.0001 arcseconds to radians
    constexpr double UNITS_2_RADIANS = ARCSECS_2_RADIANS / 1.0e4;

    // interval between basic epoch J2000.0 and current epoch (JC)
    const double centuries = (date - 51544.5) / 36525.0;

    //
    // FUNDAMENTAL ARGUMENTS in the FK5 reference system
    //

    // mean longitude of the Moon minus mean longitude of the Moon's perigee
    const double el = drange(ARCSECS_2_RADIANS * (485866.733 + (1325.0 * TURNS_2_ARCSECS + 715922.633 +
        (31.310 + 0.064 * centuries) * centuries) * centuries));

    // mean longitude of the Sun minus mean longitude of the Sun's perigee
    const double elp = drange(ARCSECS_2_RADIANS * (1287099.804 + (99.0 * TURNS_2_ARCSECS + 1292581.224 +
        (-0.577 - 0.012 * centuries) * centuries) * centuries));

    // mean longitude of the Moon minus mean longitude of the Moon's node
    const double f = drange(ARCSECS_2_RADIANS * (335778.877 + (1342.0 * TURNS_2_ARCSECS + 295263.137 +
        (-13.257 + 0.011 * centuries) * centuries) * centuries));

    // mean elongation of the Moon from the Sun
    const double d = drange(ARCSECS_2_RADIANS * (1072261.307 + (1236.0 * TURNS_2_ARCSECS + 1105601.328 +
        (-6.891 + 0.019 * centuries) * centuries) * centuries));

    // longitude of the mean ascending node of the lunar orbit on the ecliptic, measured from the mean equinox of date
    const double om = drange(ARCSECS_2_RADIANS * (450160.280 + (-5.0 * TURNS_2_ARCSECS - 482890.539 +
        (7.455 + 0.008 * centuries) * centuries) * centuries));

    // multiples of arguments
    const double el2 = el + el;
    const double el3 = el2 + el;
    const double elp2 = elp + elp;
    const double f2 = f + f;
    const double f4 = f2 + f2;
    const double d2 = d + d;
    const double d4 = d2 + d2;
    const double om2 = om + om;

    //
    // SERIES FOR THE NUTATION
    //

    double dp = 0.0;
    double de = 0.0;

    // 106
    dp = dp + std::sin(elp + d);
    // 105
    dp = dp - std::sin(f2 + d4 + om2);
    // 104
    dp = dp + std::sin(el2 + d2);
    // 103
    dp = dp - std::sin(el - f2 + d2);
    // 102
    dp = dp - std::sin(el + elp - d2 + om);
    // 101
    dp = dp - std::sin(-elp + f2 + om);
    // 100
    dp = dp - std::sin(el - f2 - d2);
    // 99
    dp = dp - std::sin(elp + d2);
    // 98
    dp = dp - std::sin(f2 - d + om2);
    // 97
    dp = dp - std::sin(-f2 + om);
    // 96
    dp = dp + std::sin(-el - elp + d2 + om);
    // 95
    dp = dp + std::sin(elp + f2 + om);
    // 94
    dp = dp - std::sin(el + f2 - d2);
    // 93
    dp = dp + std::sin(el3 + f2 - d2 + om2);
    // 92
    dp = dp + std::sin(f4 - d2 + om2);
    // 91
    dp = dp - std::sin(el + d2 + om);
    // 90
    dp = dp - std::sin(el2 + f2 + d2 + om2);
    // 89
    double angle = el2 + f2 - d2 + om;
    dp = dp + std::sin(angle);
    de = de - std::cos(angle);
    // 88
    dp = dp + std::sin(el - elp - d2);
    // 87
    dp = dp + std::sin(-el + f4 + om2);
    // 86
    angle = -el2 + f2 + d4 + om2;
    dp = dp - std::sin(angle);
    de = de + std::cos(angle);
    // 85
    angle = el + f2 + d2 + om;
    dp = dp - std::sin(angle);
    de = de + std::cos(angle);
    // 84
    angle = el + elp + f2 - d2 + om2;
    dp = dp + std::sin(angle);
    de = de - std::cos(angle);
    // 83
    dp = dp - std::sin(el2 - d4);
    // 82
    angle = -el + f2 + d4 + om2;
    dp = dp - 2.0 * std::sin(angle);
    de = de + std::cos(angle);
    // 81
    angle = -el2 + f2 + d2 + om2;
    dp = dp + std::sin(angle);
    de = de - std::cos(angle);
    // 80
    dp = dp - std::sin(el - d4);
    // 79
    angle = -el + om2;
    dp = dp + std::sin(angle);
    de = de - std::cos(angle);
    // 78
    angle = f2 + d + om2;
    dp = dp + 2.0 * std::sin(angle);
    de = de - std::cos(angle);
    // 77
    dp = dp + 2.0 * std::sin(el3);
    // 76
    angle = el + om2;
    dp = dp - 2.0 * std::sin(angle);
    de = de + std::cos(angle);
    // 75
    angle = el2 + om;
    dp = dp + 2.0 * std::sin(angle);
    de = de - std::cos(angle);
    // 74
    angle = -el + f2 - d2 + om;
    dp = dp - 2.0 * std::sin(angle);
    de = de + std::cos(angle);
    // 73
    angle = el + elp + f2 + om2;
    dp = dp + 2.0 * std::sin(angle);
    de = de - std::cos(angle);
    // 72
    angle = -elp + f2 + d2 + om2;
    dp = dp - 3.0 * std::sin(angle);
    de = de + std::cos(angle);
    // 71
    angle = el3 + f2 + om2;
    dp = dp - 3.0 * std::sin(angle);
    de = de + std::cos(angle);
    // 70
    angle = -el2 + om;
    dp = dp - 2.0 * std::sin(angle);
    de = de + std::cos(angle);
    // 69
    angle = -el - elp + f2 + d2 + om2;
    dp = dp - 3.0 * std::sin(angle);
    de = de + std::cos(angle);
    // 68
    angle = el - elp + f2 + om2;
    dp = dp - 3.0 * std::sin(angle);
    de = de + std::cos(angle);
    // 67
    dp = dp + 3.0 * std::sin(el + f2);
    // 66
    dp = dp - 3.0 * std::sin(el + elp);
    // 65
    dp = dp - 4.0 * std::sin(d);
    // 64
    dp = dp + 4.0 * std::sin(el - f2);
    // 63
    dp = dp - 4.0 * std::sin(elp - d2);
    // 62
    angle = el2 + f2 + om;
    dp = dp - 5.0 * std::sin(angle);
    de = de + 3.0 * std::cos(angle);
    // 61
    dp = dp + 5.0 * std::sin(el - elp);
    // 60
    angle = -d2 + om;
    dp = dp - 5.0 * std::sin(angle);
    de = de + 3.0 * std::cos(angle);
    // 59
    angle = el + f2 - d2 + om;
    dp = dp + 6.0 * std::sin(angle);
    de = de - 3.0 * std::cos(angle);
    // 58
    angle = f2 + d2 + om;
    dp = dp - 7.0 * std::sin(angle);
    de = de + 3.0 * std::cos(angle);
    // 57
    angle = d2 + om;
    dp = dp - 6.0 * std::sin(angle);
    de = de + 3.0 * std::cos(angle);
    // 56
    angle = el2 + f2 - d2 + om2;
    dp = dp + 6.0 * std::sin(angle);
    de = de - 3.0 * std::cos(angle);
    // 55
    dp = dp + 6.0 * std::sin(el + d2);
    // 54
    angle = el + f2 + d2 + om2;
    dp = dp - 8.0 * std::sin(angle);
    de = de + 3.0 * std::cos(angle);
    // 53
    angle = -elp + f2 + om2;
    dp = dp - 7.0 * std::sin(angle);
    de = de + 3.0 * std::cos(angle);
    // 52
    angle = elp + f2 + om2;
    dp = dp + 7.0 * std::sin(angle);
    de = de - 3.0 * std::cos(angle);
    // 51
    dp = dp - 7.0 * std::sin(el + elp - d2);
    // 50
    angle = -el + f2 + d2 + om;
    dp = dp - 10.0 * std::sin(angle);
    de = de + 5.0 * std::cos(angle);
    // 49
    angle = el - d2 + om;
    dp = dp - 13.0 * std::sin(angle);
    de = de + 7.0 * std::cos(angle);
    // 48
    angle = -el + d2 + om;
    dp = dp + 16.0 * std::sin(angle);
    de = de - 8.0 * std::cos(angle);
    // 47
    angle = -el + f2 + om;
    dp = dp + 21.0 * std::sin(angle);
    de = de - 10.0 * std::cos(angle);
    // 46
    dp = dp + 26.0 * std::sin(f2);
    de = de - std::cos(f2);
    // 45
    angle = el2 + f2 + om2;
    dp = dp - 31.0 * std::sin(angle);
    de = de + 13.0 * std::cos(angle);
    // 44
    angle = el + f2 - d2 + om2;
    dp = dp + 29.0 * std::sin(angle);
    de = de - 12.0 * std::cos(angle);
    // 43
    dp = dp + 29.0 * std::sin(el2);
    de = de - std::cos(el2);
    // 42
    angle = f2 + d2 + om2;
    dp = dp - 38.0 * std::sin(angle);
    de = de + 16.0 * std::cos(angle);
    // 41
    angle = el + f2 + om;
    dp = dp - 51.0 * std::sin(angle);
    de = de + 27.0 * std::cos(angle);
    // 40
    angle = -el + f2 + d2 + om2;
    dp = dp - 59.0 * std::sin(angle);
    de = de + 26.0 * std::cos(angle);
    // 39
    angle = -el + om;
    dp = dp + (-58.0 - 0.1 * centuries) * std::sin(angle);
    de = de + 32.0 * std::cos(angle);
    // 38
    angle = el + om;
    dp = dp + (63.0 + 0.1 * centuries) * std::sin(angle);
    de = de - 33.0 * std::cos(angle);
    // 37
    dp = dp + 63.0 * std::sin(d2);
    de = de - 2.0 * std::cos(d2);
    // 36
    angle = -el + f2 + om2;
    dp = dp + 123.0 * std::sin(angle);
    de = de - 53.0 * std::cos(angle);
    // 35
    angle = el - d2;
    dp = dp - 158.0 * std::sin(angle);
    de = de - std::cos(angle);
    // 34
    angle = el + f2 + om2;
    dp = dp - 301.0 * std::sin(angle);
    de = de + (129.0 - 0.1 * centuries) * std::cos(angle);
    // 33
    angle = f2 + om;
    dp = dp + (-386.0 - 0.4 * centuries) * std::sin(angle);
    de = de + 200.0 * std::cos(angle);
    // 32
    dp = dp + (712.0 + 0.1 * centuries) * std::sin(el);
    de = de - 7.0 * std::cos(el);
    // 31
    angle = f2 + om2;
    dp = dp + (-2274.0 - 0.2 * centuries) * std::sin(angle);
    de = de + (977.0 - 0.5 * centuries) * std::cos(angle);
    // 30
    dp = dp - std::sin(elp + f2 - d2);
    // 29
    dp = dp + std::sin(-el + d + om);
    // 28
    dp = dp + std::sin(elp + om2);
    // 27
    dp = dp - std::sin(elp - f2 + d2);
    // 26
    dp = dp + std::sin(-f2 + d2 + om);
    // 25
    dp = dp + std::sin(el2 + elp - d2);
    // 24
    dp = dp - 4.0 * std::sin(el - d);
    // 23
    angle = elp + f2 - d2 + om;
    dp = dp + 4.0 * std::sin(angle);
    de = de - 2.0 * std::cos(angle);
    // 22
    angle = el2 - d2 + om;
    dp = dp + 4.0 * std::sin(angle);
    de = de - 2.0 * std::cos(angle);
    // 21
    angle = -elp + f2 - d2 + om;
    dp = dp - 5.0 * std::sin(angle);
    de = de + 3.0 * std::cos(angle);
    // 20
    angle = -el2 + d2 + om;
    dp = dp - 6.0 * std::sin(angle);
    de = de + 3.0 * std::cos(angle);
    // 19
    angle = -elp + om;
    dp = dp - 12.0 * std::sin(angle);
    de = de + 6.0 * std::cos(angle);
    // 18
    angle = elp2 + f2 - d2 + om2;
    dp = dp + (-16.0 + 0.1 * centuries) * std::sin(angle);
    de = de + 7.0 * std::cos(angle);
    // 17
    angle = elp + om;
    dp = dp - 15.0 * std::sin(angle);
    de = de + 9.0 * std::cos(angle);
    // 16
    dp = dp + (17.0 - 0.1 * centuries) * std::sin(elp2);
    // 15
    dp = dp - 22.0 * std::sin(f2 - d2);
    // 14
    angle = el2 - d2;
    dp = dp + 48.0 * std::sin(angle);
    de = de + std::cos(angle);
    // 13
    angle = f2 - d2 + om;
    dp = dp + (129.0 + 0.1 * centuries) * std::sin(angle);
    de = de - 70.0 * std::cos(angle);
    // 12
    angle = -elp + f2 - d2 + om2;
    dp = dp + (217.0 - 0.5 * centuries) * std::sin(angle);
    de = de + (-95.0 + 0.3 * centuries) * std::cos(angle);
    // 11
    angle = elp + f2 - d2 + om2;
    dp = dp + (-517.0 + 1.2 * centuries) * std::sin(angle);
    de = de + (224.0 - 0.6 * centuries) * std::cos(angle);
    // 10
    dp = dp + (1426.0 - 3.4 * centuries) * std::sin(elp);
    de = de + (54.0 - 0.1 * centuries) * std::cos(elp);
    // 9
    angle = f2 - d2 + om2;
    dp = dp + (-13187.0 - 1.6 * centuries) * std::sin(angle);
    de = de + (5736.0 - 3.1 * centuries) * std::cos(angle);
    // 8
    dp = dp + std::sin(el2 - f2 + om);
    // 7
    angle = -elp2 + f2 - d2 + om;
    dp = dp - 2.0 * std::sin(angle);
    de = de + 1.0 * std::cos(angle);
    // 6
    dp = dp - 3.0 * std::sin(el - elp - d);
    // 5
    angle = -el2 + f2 + om2;
    dp = dp - 3.0 * std::sin(angle);
    de = de + 1.0 * std::cos(angle);
    // 4
    dp = dp + 11.0 * std::sin(el2 - f2);
    // 3
    angle = -el2 + f2 + om;
    dp = dp + 46.0 * std::sin(angle);
    de = de - 24.0 * std::cos(angle);
    // 2
    dp = dp + (2062.0 + 0.2 * centuries) * std::sin(om2);
    de = de + (-895.0 + 0.5 * centuries) * std::cos(om2);
    // 1
    dp = dp + (-171996.0 - 174.2 * centuries) * std::sin(om);
    de = de + (92025.0 + 8.9 * centuries) * std::cos(om);

    // convert results to radians
    psi = dp * UNITS_2_RADIANS;
    eps = de * UNITS_2_RADIANS;

    // mean obliquity
    eps0 = ARCSECS_2_RADIANS * (84381.448 + (-46.8150 + (-0.00059 + 0.001813 * centuries) * centuries) * centuries);
}

}
