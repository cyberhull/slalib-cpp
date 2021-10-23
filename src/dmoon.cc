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
 * Calculates approximate geocentric position and velocity of the Moon (double precision).
 *
 * This function is a full implementation of the algorithm published by Meeus (see reference).
 *
 * Meeus quotes accuracies of 10 arcsec in longitude, 3 arcsec in latitude and 0.2 arcsec in HP (equivalent to about
 * 20 km in distance). Comparison with JPL DE200 over the interval 1960-2025 gives RMS errors of 3.7 arcsec and 83
 * mas/hour in longitude, 2.3 arcsec and 48 mas/hour in latitude, 11 km and 81 mm/s in distance. The maximum errors
 * over the same interval are 18 arcsec and 0.50 arcsec/hour in longitude, 11 arcsec and 0.24 arcsec/hour in
 * latitude, 40 km and 0.29 m/s in distance.
 *
 * The original algorithm is expressed in terms of the obsolete timescale Ephemeris Time. Either TDB or TT can be
 * used, but not UT without incurring significant errors (30 arcsec at the present time) due to the Moon's 0.5
 * arcsec/sec movement.
 *
 * The algorithm is based on pre IAU 1976 standards. However, the result has been moved onto the new (FK5) equinox, an
 * adjustment which is in any case much smaller than the intrinsic accuracy of the procedure.
 *
 * Velocity is obtained by a complete analytical differentiation of the Meeus model.
 *
 * Reference:
 *    Meeus, l'Astronomie, June 1984, p. 348.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param date TDB (Barycentric Dynamical Time; loosely ET) as a Modified Julian Date (JD-2400000.5).
 * @param pv Return value: Moon {x,y,z},{xdot,ydot,zdot}, mean equator and equinox of date (AU, AU/s).
 */
void dmoon(double date, VectorPV<double>& pv) {
    // degrees, arcseconds and seconds of time to radians
    constexpr double DEGREES_2_RADIANS = 0.0174532925199432957692369;
    constexpr double ARCSECONDS_2_RADIANS = 4.848136811095359935899141e-6;
    constexpr double SECONDS_2_RADIANS = 7.272205216643039903848712e-5;

    // seconds per Julian century (86400*36525)
    constexpr double SECONDS_PER_JCENTURY = 3155760000.0;

    // Julian epoch of JEPOCH_B1950
    constexpr double JEPOCH_B1950 = 1949.9997904423;

    // earth equatorial radius in AU ( = 6378.137 / 149597870 )
    constexpr double EARTH_RADIUS_AU = 4.2635212653763e-5;

    /*
     * Coefficients for fundamental arguments
     *
     *   at J1900:  t**0, t**1, t**2, t**3
     *   at epoch:  t**0, t**1
     *
     * Units are degrees for position and Julian centuries for time.
     */

    // Moon's mean longitude
    constexpr double ELP0 = 270.434164;
    constexpr double ELP1 = 481267.8831;
    constexpr double ELP2 = -0.001133;
    constexpr double ELP3 = 0.0000019;
    double elp, delp;

    // Sun's mean anomaly
    constexpr double EM0 = 358.475833;
    constexpr double EM1 = 35999.0498;
    constexpr double EM2 = -0.000150;
    constexpr double EM3 = -0.0000033;
    double em, dem;

    // Moon's mean anomaly
    constexpr double EMP0 = 296.104608;
    constexpr double EMP1 = 477198.8491;
    constexpr double EMP2 = 0.009192;
    constexpr double EMP3 = 0.0000144;
    double emp, demp;

    // Moon's mean elongation
    constexpr double D0 = 350.737486;
    constexpr double D1 = 445267.1142;
    constexpr double D2 = -0.001436;
    constexpr double D3 = 0.0000019;
    double d, dd;

    // mean distance of the Moon from its ascending node
    constexpr double F0 = 11.250889;
    constexpr double F1 = 483202.0251;
    constexpr double F2 = -0.003211;
    constexpr double F3 = -0.0000003;
    double f, df;

    // longitude of the Moon's ascending node
    constexpr double OM0 = 259.183275;
    constexpr double OM1 = -1934.1420;
    constexpr double OM2 = 0.002078;
    constexpr double OM3 = 0.0000022;
    double om, dom;

    // coefficients for (dimensionless) e factor
    constexpr double E1 = -0.002495;
    constexpr double E2 = -0.00000752;
    double e, de, esq, desq;

    // coefficients for periodic variations etc.
    constexpr double PAC = 0.000233;
    constexpr double PA0 = 51.2;
    constexpr double PA1 = 20.2;
    constexpr double PBC = -0.001778;
    constexpr double PCC = 0.000817;
    constexpr double PDC = 0.002011;
    constexpr double PEC = 0.003964;
    constexpr double PE0 = 346.560;
    constexpr double PE1 = 132.870;
    constexpr double PE2 = -0.0091731;
    constexpr double PFC = 0.001964;
    constexpr double PGC = 0.002541;
    constexpr double PHC = 0.001964;
    constexpr double PIC = -0.024691;
    constexpr double PJC = -0.004328;
    constexpr double PJ0 = 275.05;
    constexpr double PJ1 = -2.30;
    constexpr double CW1 = 0.0004664;
    constexpr double CW2 = 0.0000754;

    /*
     *  Coefficients for Moon position
     *
     *   Tx(n)       = coefficient of L, b or p term (deg)
     *   ITx(n,1-5)  = coefficients of M, M', d, f, e**n in argument
     */

    // longitude
    constexpr int NL = 50;
    static const double TL[NL] = {
        +6.288750, +1.274018, +0.658309, +0.213616, -0.185596,
        -0.114336, +0.058793, +0.057212, +0.053320, +0.045874,
        +0.041024, -0.034718, -0.030465, +0.015326, -0.012528,
        -0.010980, +0.010674, +0.010034, +0.008548, -0.007910,
        -0.006783, +0.005162, +0.005000, +0.004049, +0.003996,
        +0.003862, +0.003665, +0.002695, +0.002602, +0.002396,
        -0.002349, +0.002249, -0.002125, -0.002079, +0.002059,
        -0.001773, -0.001595, +0.001220, -0.001110, +0.000892,
        -0.000811, +0.000761, +0.000717, +0.000704, +0.000693,
        +0.000598, +0.000550, +0.000538, +0.000521, +0.000486
    };
    static const int8_t ITL[NL][5] = {
    //    M   M'  d   f  n
        {+0, +1, +0, +0, 0},
        {+0, -1, +2, +0, 0},
        {+0, +0, +2, +0, 0},
        {+0, +2, +0, +0, 0},
        {+1, +0, +0, +0, 1},
        {+0, +0, +0, +2, 0},
        {+0, -2, +2, +0, 0},
        {-1, -1, +2, +0, 1},
        {+0, +1, +2, +0, 0},
        {-1, +0, +2, +0, 1},
        {-1, +1, +0, +0, 1},
        {+0, +0, +1, +0, 0},
        {+1, +1, +0, +0, 1},
        {+0, +0, +2, -2, 0},
        {+0, +1, +0, +2, 0},
        {+0, -1, +0, +2, 0},
        {+0, -1, +4, +0, 0},
        {+0, +3, +0, +0, 0},
        {+0, -2, +4, +0, 0},
        {+1, -1, +2, +0, 1},
        {+1, +0, +2, +0, 1},
        {+0, +1, -1, +0, 0},
        {+1, +0, +1, +0, 1},
        {-1, +1, +2, +0, 1},
        {+0, +2, +2, +0, 0},
        {+0, +0, +4, +0, 0},
        {+0, -3, +2, +0, 0},
        {-1, +2, +0, +0, 1},
        {+0, +1, -2, -2, 0},
        {-1, -2, +2, +0, 1},
        {+0, +1, +1, +0, 0},
        {-2, +0, +2, +0, 2},
        {+1, +2, +0, +0, 1},
        {+2, +0, +0, +0, 2},
        {-2, -1, +2, +0, 2},
        {+0, +1, +2, -2, 0},
        {+0, +0, +2, +2, 0},
        {-1, -1, +4, +0, 1},
        {+0, +2, +0, +2, 0},
        {+0, +1, -3, +0, 0},
        {+1, +1, +2, +0, 1},
        {-1, -2, +4, +0, 1},
        {-2, +1, +0, +0, 2},
        {-2, +1, -2, +0, 2},
        {+1, -2, +2, +0, 1},
        {-1, +0, +2, -2, 1},
        {+0, +1, +4, +0, 0},
        {+0, +4, +0, +0, 0},
        {-1, +0, +4, +0, 1},
        {+0, +2, -1, +0, 0}
    };

    // latitude
    constexpr int NB = 45;
    static const double TB[NB] = {
        +5.128189, +0.280606, +0.277693, +0.173238, +0.055413,
        +0.046272, +0.032573, +0.017198, +0.009267, +0.008823,
        +0.008247, +0.004323, +0.004200, +0.003372, +0.002472,
        +0.002222, +0.002072, +0.001877, +0.001828, -0.001803,
        -0.001750, +0.001570, -0.001487, -0.001481, +0.001417,
        +0.001350, +0.001330, +0.001106, +0.001020, +0.000833,
        +0.000781, +0.000670, +0.000606, +0.000597, +0.000492,
        +0.000450, +0.000439, +0.000423, +0.000422, -0.000367,
        -0.000353, +0.000331, +0.000317, +0.000306, -0.000283
    };
    static const int8_t ITB[NB][5] = {
    //    M   M'  d   f  n
        {+0, +0, +0, +1, 0},
        {+0, +1, +0, +1, 0},
        {+0, +1, +0, -1, 0},
        {+0, +0, +2, -1, 0},
        {+0, -1, +2, +1, 0},
        {+0, -1, +2, -1, 0},
        {+0, +0, +2, +1, 0},
        {+0, +2, +0, +1, 0},
        {+0, +1, +2, -1, 0},
        {+0, +2, +0, -1, 0},
        {-1, +0, +2, -1, 1},
        {+0, -2, +2, -1, 0},
        {+0, +1, +2, +1, 0},
        {-1, +0, -2, +1, 1},
        {-1, -1, +2, +1, 1},
        {-1, +0, +2, +1, 1},
        {-1, -1, +2, -1, 1},
        {-1, +1, +0, +1, 1},
        {+0, -1, +4, -1, 0},
        {+1, +0, +0, +1, 1},
        {+0, +0, +0, +3, 0},
        {-1, +1, +0, -1, 1},
        {+0, +0, +1, +1, 0},
        {+1, +1, +0, +1, 1},
        {-1, -1, +0, +1, 1},
        {-1, +0, +0, +1, 1},
        {+0, +0, -1, +1, 0},
        {+0, +3, +0, +1, 0},
        {+0, +0, +4, -1, 0},
        {+0, -1, +4, +1, 0},
        {+0, +1, +0, -3, 0},
        {+0, -2, +4, +1, 0},
        {+0, +0, +2, -3, 0},
        {+0, +2, +2, -1, 0},
        {-1, +1, +2, -1, 1},
        {+0, +2, -2, -1, 0},
        {+0, +3, +0, -1, 0},
        {+0, +2, +2, +1, 0},
        {+0, -3, +2, -1, 0},
        {+1, -1, +2, +1, 1},
        {+1, +0, +2, +1, 1},
        {+0, +0, +4, +1, 0},
        {-1, +1, +2, +1, 1},
        {-2, +0, +2, -1, 2},
        {+0, +1, +0, +3, 0}
    };

    // parallax
    constexpr int NP = 31;
    double TP[NP]{
        +0.950724, +0.051818, +0.009531, +0.007843, +0.002824,
        +0.000857, +0.000533, +0.000401, +0.000320, -0.000271,
        -0.000264, -0.000198, +0.000173, +0.000167, -0.000111,
        +0.000103, -0.000084, -0.000083, +0.000079, +0.000072,
        +0.000064, -0.000063, +0.000041, +0.000035, -0.000033,
        -0.000030, -0.000029, -0.000029, +0.000026, -0.000023,
        +0.000019
    };
    static const int8_t ITP[NP][5] = {
    //    M   M'  d   f  n
        {+0, +0, +0, +0, 0},
        {+0, +1, +0, +0, 0},
        {+0, -1, +2, +0, 0},
        {+0, +0, +2, +0, 0},
        {+0, +2, +0, +0, 0},
        {+0, +1, +2, +0, 0},
        {-1, +0, +2, +0, 1},
        {-1, -1, +2, +0, 1},
        {-1, +1, +0, +0, 1},
        {+0, +0, +1, +0, 0},
        {+1, +1, +0, +0, 1},
        {+0, -1, +0, +2, 0},
        {+0, +3, +0, +0, 0},
        {+0, -1, +4, +0, 0},
        {+1, +0, +0, +0, 1},
        {+0, -2, +4, +0, 0},
        {+0, +2, -2, +0, 0},
        {+1, +0, +2, +0, 1},
        {+0, +2, +2, +0, 0},
        {+0, +0, +4, +0, 0},
        {-1, +1, +2, +0, 1},
        {+1, -1, +2, +0, 1},
        {+1, +0, +1, +0, 1},
        {-1, +2, +0, +0, 1},
        {+0, +3, -2, +0, 0},
        {+0, +1, +1, +0, 0},
        {+0, +0, -2, +2, 0},
        {+1, +2, +0, +0, 1},
        {-2, +0, +2, +0, 2},
        {+0, +1, -2, +2, 0},
        {-1, -1, +4, +0, 1}
    };

    // centuries since J1900
    double t = (date - 15019.5) / 36525.0;

    /*
     * Fundamental arguments (radians) and derivatives (radians per
     * Julian century) for the current epoch.
     */

    // Moon's mean longitude
    elp = DEGREES_2_RADIANS * std::fmod(ELP0 + (ELP1 + (ELP2 + ELP3 * t) * t) * t, 360.0);
    delp = DEGREES_2_RADIANS * (ELP1 + (2.0 * ELP2 + 3.0 * ELP3 * t) * t);

    // Sun's mean anomaly
    em = DEGREES_2_RADIANS * std::fmod(EM0 + (EM1 + (EM2 + EM3 * t) * t) * t, 360.0);
    dem = DEGREES_2_RADIANS * (EM1 + (2.0 * EM2 + 3.0 * EM3 * t) * t);

    // Moon's mean anomaly
    emp = DEGREES_2_RADIANS * std::fmod(EMP0 + (EMP1 + (EMP2 + EMP3 * t) * t) * t, 360.0);
    demp = DEGREES_2_RADIANS * (EMP1 + (2.0 * EMP2 + 3.0 * EMP3 * t) * t);

    // Moon's mean elongation
    d = DEGREES_2_RADIANS * std::fmod(D0 + (D1 + (D2 + D3 * t) * t) * t, 360.0);
    dd = DEGREES_2_RADIANS * (D1 + (2.0 * D2 + 3.0 * D3 * t) * t);

    // mean distance of the Moon from its ascending node
    f = DEGREES_2_RADIANS * std::fmod(F0 + (F1 + (F2 + F3 * t) * t) * t, 360.0);
    df = DEGREES_2_RADIANS * (F1 + (2.0 * F2 + 3.0 * F3 * t) * t);

    // longitude of the Moon's ascending node
    om = DEGREES_2_RADIANS * std::fmod(OM0 + (OM1 + (OM2 + OM3 * t) * t) * t, 360.0);
    dom = DEGREES_2_RADIANS * (OM1 + (2.0 * OM2 + 3.0 * OM3 * t) * t);
    const double sin_om = std::sin(om);
    const double cos_om = std::cos(om);
    const double dom_cos_om = dom * cos_om;

    // add the periodic variations
    double theta = DEGREES_2_RADIANS * (PA0 + PA1 * t);
    const double wa = std::sin(theta);
    const double dwa = DEGREES_2_RADIANS * PA1 * std::cos(theta);
    theta = DEGREES_2_RADIANS * (PE0 + (PE1 + PE2 * t) * t);
    const double wb = PEC * std::sin(theta);
    const double dwb = DEGREES_2_RADIANS * PEC * (PE1 + 2.0 * PE2 * t) * std::cos(theta);
    elp = elp + DEGREES_2_RADIANS * (PAC * wa + wb + PFC * sin_om);
    delp = delp + DEGREES_2_RADIANS * (PAC * dwa + dwb + PFC * dom_cos_om);
    em = em + DEGREES_2_RADIANS * PBC * wa;
    dem = dem + DEGREES_2_RADIANS * PBC * dwa;
    emp = emp + DEGREES_2_RADIANS * (PCC * wa + wb + PGC * sin_om);
    demp = demp + DEGREES_2_RADIANS * (PCC * dwa + dwb + PGC * dom_cos_om);
    d = d + DEGREES_2_RADIANS * (PDC * wa + wb + PHC * sin_om);
    dd = dd + DEGREES_2_RADIANS * (PDC * dwa + dwb + PHC * dom_cos_om);
    const double wom = om + DEGREES_2_RADIANS * (PJ0 + PJ1 * t);
    const double dwom = dom + DEGREES_2_RADIANS * PJ1;
    const double sin_wom = std::sin(wom);
    const double cos_wom = std::cos(wom);
    f = f + DEGREES_2_RADIANS * (wb + PIC * sin_om + PJC * sin_wom);
    df = df + DEGREES_2_RADIANS * (dwb + PIC * dom_cos_om + PJC * dwom * cos_wom);

    // e-factor, and square
    e = 1.0 + (E1 + E2 * t) * t;
    de = E1 + 2.0 * E2 * t;
    esq = e * e;
    desq = 2.0 * e * de;

    /*
     *  Series expansions.
     */
    int n, i;

    // longitude
    double v = 0.0;
    double dv = 0.0;
    double coeff, emn, empn, dn, fn, en, den, dtheta, ftheta;
    for (n = NL - 1; n >= 0; n--) {
        coeff = TL[n];
        emn = (double) ITL[n][0];
        empn = (double) ITL[n][1];
        dn = (double) ITL[n][2];
        fn = (double) ITL[n][3];
        i = ITL[n][4];
        if (i == 0) {
            en = 1.0;
            den = 0.0;
        } else if (i == 1) {
            en = e;
            den = de;
        } else {
            en = esq;
            den = desq;
        }
        theta = emn * em + empn * emp + dn * d + fn * f;
        dtheta = emn * dem + empn * demp + dn * dd + fn * df;
        ftheta = std::sin(theta);
        v += coeff * ftheta * en;
        dv += coeff * (std::cos(theta) * dtheta * en + ftheta * den);
    }
    const double el = elp + DEGREES_2_RADIANS * v;
    const double del = (delp + DEGREES_2_RADIANS * dv) / SECONDS_PER_JCENTURY;

    // latitude
    v = 0.0;
    dv = 0.0;
    for (n = NB - 1; n >= 0; n--) {
        coeff = TB[n];
        emn = (double) ITB[n][0];
        empn = (double) ITB[n][1];
        dn = (double) ITB[n][2];
        fn = (double) ITB[n][3];
        i = ITB[n][4];
        if (i == 0) {
            en = 1.0;
            den = 0.0;
        } else if (i == 1) {
            en = e;
            den = de;
        } else {
            en = esq;
            den = desq;
        }
        theta = emn * em + empn * emp + dn * d + fn * f;
        dtheta = emn * dem + empn * demp + dn * dd + fn * df;
        ftheta = std::sin(theta);
        v += coeff * ftheta * en;
        dv += coeff * (std::cos(theta) * dtheta * en + ftheta * den);
    }
    const double bf = 1.0 - CW1 * cos_om - CW2 * cos_wom;
    const double dbf = CW1 * dom * sin_om + CW2 * dwom * sin_wom;
    const double b = DEGREES_2_RADIANS * v * bf;
    const double db = DEGREES_2_RADIANS * (dv * bf + v * dbf) / SECONDS_PER_JCENTURY;

    // parallax
    v = 0.0;
    dv = 0.0;
    for (n = NP - 1; n >= 0; n--) {
        coeff = TP[n];
        emn = (double) ITP[n][0];
        empn = (double) ITP[n][1];
        dn = (double) ITP[n][2];
        fn = (double) ITP[n][3];
        i = ITP[n][4];
        if (i == 0) {
            en = 1.0;
            den = 0.0;
        } else if (i == 1) {
            en = e;
            den = de;
        } else {
            en = esq;
            den = desq;
        }
        theta = emn * em + empn * emp + dn * d + fn * f;
        dtheta = emn * dem + empn * demp + dn * dd + fn * df;
        ftheta = std::cos(theta);
        v += coeff * ftheta * en;
        dv += coeff * (-std::sin(theta) * dtheta * en + ftheta * den);
    }
    const double p = DEGREES_2_RADIANS * v;
    const double dp = DEGREES_2_RADIANS * dv / SECONDS_PER_JCENTURY;

    /*
     * Transformation into final form.
     */

    // parallax to distance (AU, AU/sec)
    const double sin_p = std::sin(p);
    const double r = EARTH_RADIUS_AU / sin_p;
    const double dr = -r * dp * std::cos(p) / sin_p;

    // longitude, latitude to x,y,z (AU)
    const double sel = std::sin(el);
    const double cel = std::cos(el);
    const double sb = std::sin(b);
    const double cb = std::cos(b);
    const double rcb = r * cb;
    const double rbd = r * db;
    const double w = rbd * sb - cb * dr;
    const double x = rcb * cel;
    const double y = rcb * sel;
    const double z = r * sb;
    const double xd = -y * del - w * cel;
    const double yd = x * del - w * sel;
    const double zd = rbd * cb + sb * dr;

    // Julian centuries since J2000
    t = (date - 51544.5) / 36525.0;

    // Fricke equinox correction
    const double epj = 2000.0 + t * 100.0;
    const double eqcor = SECONDS_2_RADIANS * (0.035 + 0.00085 * (epj - JEPOCH_B1950));

    // mean obliquity (IAU 1976)
    const double eps = ARCSECONDS_2_RADIANS * (84381.448 + (-46.8150 + (-0.00059 + 0.001813 * t) * t) * t);

    // to the equatorial system, mean of date, FK5 system
    const double sin_eps = std::sin(eps);
    const double cos_eps = std::cos(eps);
    const double es = eqcor * sin_eps;
    const double ec = eqcor * cos_eps;
    pv.set_x(x - ec * y + es * z);
    pv.set_y(eqcor * x + y * cos_eps - z * sin_eps);
    pv.set_z(y * sin_eps + z * cos_eps);
    pv.set_dx(xd - ec * yd + es * zd);
    pv.set_dy(eqcor * xd + yd * cos_eps - zd * sin_eps);
    pv.set_dz(yd * sin_eps + zd * cos_eps);
}

}
