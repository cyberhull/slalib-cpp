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
#include <cstdint>

namespace sla {

/**
 * Calculates approximate geocentric position and velocity of the Moon (single precision).
 *
 * The date and time is TDB (Barycentric Dynamical Time; loosely ET) in a Julian calendar, which has been aligned to
 * the ordinary Gregorian calendar for the interval 1900 March 1 to 2100 February 28. The `year` and `day` can be
 * obtained by calling sla::calyd() or sla::clyd() functions.
 *
 * The returned position is accurate to better than 0.5 arcminute in direction and 1000 km in distance. The
 * returned velocity is accurate to better than 0.5"/hour in direction and 4 m/s in distance. (RMS figures with
 * respect to JPL DE200 for the interval 1960-2025 are 14 arcsec and 0.2 arcsec/hour in longitude, 9 arcsec and 0.2
 * arcsec/hour in latitude, 350 km and 2 m/s in distance.) Note that the distance accuracy is comparatively poor
 * because this function is principally intended for computing topocentric direction.
 *
 * This function is only a partial implementation of the original Meeus algorithm (reference below), which offers 4
 * times the accuracy in direction and 30 times the accuracy in distance when fully implemented (as it is in the
 * sla::dmoon() function).
 *
 * Reference:
 *   Meeus, l'Astronomie, June 1984, p348.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param year Year.
 * @param day Day in year (1 = Jan 1-st).
 * @param fraction Fraction of day.
 * @param pv Return value: Moon position and velocity vector: Moon center relative to Earth center, mean equator and
 *   equinox of date; position part is in AU; velocity part is in AU/sec.
 */
void moon(int year, int day, float fraction, VectorPV<float>& pv) {
    // degrees to radians
    constexpr float DEGREES_2_RADIANS = 1.745329252e-2f;

    // rate conversion factor: DEGREES_2_RADIANS**2/(86400*365.25)
    constexpr float RATE_CONV_FACTOR = 9.652743551e-12f;

    // Earth radius in AU: 6378.137/149597870
    constexpr float EARTH_RADIUS_AU = 4.2635212653763e-5f;

    /*
     * Coefficients for fundamental arguments:
     *
     * fixed term (deg), term in t (deg & whole revs + fraction per year).
     */

    // Moon's mean longitude
    constexpr float ELP0 = 270.434164f,
        ELP1 = 4812.678831f,
        ELP1I = 4680.0f,
        ELP1F = 132.678831f;

    // Sun's mean anomaly
    constexpr float EM0 = 358.475833f,
        EM1 = 359.990498f,
        EM1F = 359.990498f;

    // Moon's mean anomaly
    constexpr float EMP0 = 296.104608f,
        EMP1 = 4771.988491f,
        EMP1I = 4680.0f,
        EMP1F = 91.988491f;

    // Moon's mean elongation
    constexpr float D0 = 350.737486f,
        D1 = 4452.671142f,
        D1I = 4320.0f,
        D1F = 132.671142f;

    // mean distance of the Moon from its ascending node
    constexpr float F0 = 11.250889f,
        F1 = 4832.020251f,
        F1I = 4680.0f,
        F1F = 152.020251;

    /*
     * Coefficients for Moon position:
     *
     * t[n]        = coefficient of term (deg.)
     * IT[n][0..3] = coefficients of M, M', d, f in argument
     */

    // longitude
    static const float TL[39] = {
        +6.288750f, +1.274018f, +0.658309f, +0.213616f, -0.185596f,
        -0.114336f, +0.058793f, +0.057212f, +0.053320f, +0.045874f,
        +0.041024f, -0.034718f, -0.030465f, +0.015326f, -0.012528f,
        -0.010980f, +0.010674f, +0.010034f, +0.008548f, -0.007910f,
        -0.006783f, +0.005162f, +0.005000f, +0.004049f, +0.003996f,
        +0.003862f, +0.003665f, +0.002695f, +0.002602f, +0.002396f,
        -0.002349f, +0.002249f, -0.002125f, -0.002079f, +0.002059f,
        -0.001773f, -0.001595f, +0.001220f, -0.001110f
    };
    static const int8_t ITL[39][4] = {
        //   M   M'  d   f
        {0,  +1, 0,  0},
        {0,  -1, +2, 0},
        {0,  0,  +2, 0},
        {0,  +2, 0,  0},
        {+1, 0,  0,  0},
        {0,  0,  0,  +2},
        {0,  -2, +2, 0},
        {-1, -1, +2, 0},
        {0,  +1, +2, 0},
        {-1, 0,  +2, 0},
        {-1, +1, 0,  0},
        {0,  0,  +1, 0},
        {+1, +1, 0,  0},
        {0,  0,  +2, -2},
        {0,  +1, 0,  +2},
        {0,  -1, 0,  +2},
        {0,  -1, +4, 0},
        {0,  +3, 0,  0},
        {0,  -2, +4, 0},
        {+1, -1, +2, 0},
        {+1, 0,  +2, 0},
        {0,  +1, -1, 0},
        {+1, 0,  +1, 0},
        {-1, +1, +2, 0},
        {0,  +2, +2, 0},
        {0,  0,  +4, 0},
        {0,  -3, +2, 0},
        {-1, +2, 0,  0},
        {0,  +1, -2, -2},
        {-1, -2, +2, 0},
        {0,  +1, +1, 0},
        {-2, 0,  +2, 0},
        {+1, +2, 0,  0},
        {+2, 0,  0,  0},
        {-2, -1, +2, 0},
        {0,  +1, +2, -2},
        {0,  0,  +2, +2},
        {-1, -1, +4, 0},
        {0,  +2, 0,  +2}
    };

    // latitude
    static const float TB[29] = {
        +5.128189f, +0.280606f, +0.277693f, +0.173238f, +0.055413f,
        +0.046272f, +0.032573f, +0.017198f, +0.009267f, +0.008823f,
        +0.008247f, +0.004323f, +0.004200f, +0.003372f, +0.002472f,
        +0.002222f, +0.002072f, +0.001877f, +0.001828f, -0.001803f,
        -0.001750f, +0.001570f, -0.001487f, -0.001481f, +0.001417f,
        +0.001350f, +0.001330f, +0.001106f, +0.001020f
    };
    static const int8_t ITB[29][4] = {
        //   M   M'  d   f
        {0,  0,  0,  +1},
        {0,  +1, 0,  +1},
        {0,  +1, 0,  -1},
        {0,  0,  +2, -1},
        {0,  -1, +2, +1},
        {0,  -1, +2, -1},
        {0,  0,  +2, +1},
        {0,  +2, 0,  +1},
        {0,  +1, +2, -1},
        {0,  +2, 0,  -1},
        {-1, 0,  +2, -1},
        {0,  -2, +2, -1},
        {0,  +1, +2, +1},
        {-1, 0,  -2, +1},
        {-1, -1, +2, +1},
        {-1, 0,  +2, +1},
        {-1, -1, +2, -1},
        {-1, +1, 0,  +1},
        {0,  -1, +4, -1},
        {+1, 0,  0,  +1},
        {0,  0,  0,  +3},
        {-1, +1, 0,  -1},
        {0,  0,  +1, +1},
        {+1, +1, 0,  +1},
        {-1, -1, 0,  +1},
        {-1, 0,  0,  +1},
        {0,  0,  -1, +1},
        {0,  +3, 0,  +1},
        {0,  0,  +4, -1}
    };

    // parallax
    static const float TP[4] = {
        +0.051818f, +0.009531f, +0.007843f, +0.002824f
    };
    static const int8_t ITP[4][4] = {
        //   M   M'  d   f
        {0, +1, 0,  0},
        {0, -1, +2, 0},
        {0, 0,  +2, 0},
        {0, +2, 0,  0}
    };

    // whole years & fraction of year, and years since J1900.0
    const float yi = float(year - 1900);
    const int iy4 = ((year % 4) + 4) % 4;
    const float yf = (float(4 * (day - 1 / (iy4 + 1)) - iy4 - 2) + 4.0f * fraction) / 1461.0f;
    const float t = yi + yf;

    // Moon's mean longitude
    const float elp = DEGREES_2_RADIANS * std::fmod(ELP0 + ELP1I * yf + ELP1F * t, 360.0f);

    // Sun's mean anomaly
    const float em = DEGREES_2_RADIANS * std::fmod(EM0 + EM1F * t, 360.0f);

    // Moon's mean anomaly
    const float emp = DEGREES_2_RADIANS * std::fmod(EMP0 + EMP1I * yf + EMP1F * t, 360.0f);

    // Moon's mean elongation
    const float d = DEGREES_2_RADIANS * std::fmod(D0 + D1I * yf + D1F * t, 360.0f);

    // mean distance of the moon from its ascending node
    const float f = DEGREES_2_RADIANS * std::fmod(F0 + F1I * yf + F1F * t, 360.0f);

    // longitude
    float el = 0.0f;
    float eld = 0.0f;
    float coeff, cem, cemp, cd, cf, theta, thetad;
    int i;
    for (i = 38; i >= 0; i--) {
        coeff = TL[i];
        cem = (float) ITL[i][0];
        cemp = (float) ITL[i][1];
        cd = (float) ITL[i][2];
        cf = (float) ITL[i][3];
        theta = cem * em + cemp * emp + cd * d + cf * f;
        thetad = cem * EM1 + cemp * EMP1 + cd * D1 + cf * F1;
        el += coeff * std::sin(theta);
        eld += coeff * std::cos(theta) * thetad;
    }
    el = el * DEGREES_2_RADIANS + elp;
    eld = RATE_CONV_FACTOR * (eld + ELP1 / DEGREES_2_RADIANS);

    // latitude
    float b = 0.0f;
    float bd = 0.0f;
    for (i = 28; i >= 0; i--) {
        coeff = TB[i];
        cem = (float) ITB[i][0];
        cemp = (float) ITB[i][1];
        cd = (float) ITB[i][2];
        cf = (float) ITB[i][3];
        theta = cem * em + cemp * emp + cd * d + cf * f;
        thetad = cem * EM1 + cemp * EMP1 + cd * D1 + cf * F1;
        b += coeff * std::sin(theta);
        bd += coeff * std::cos(theta) * thetad;
    }
    b *= DEGREES_2_RADIANS;
    bd *= RATE_CONV_FACTOR;

    // parallax
    float p = 0.0f;
    float pd = 0.0f;
    for (i = 3; i >= 0; i--) {
        coeff = TP[i];
        cem = (float) ITP[i][0];
        cemp = (float) ITP[i][1];
        cd = (float) ITP[i][2];
        cf = (float) ITP[i][3];
        theta = cem * em + cemp * emp + cd * d + cf * f;
        thetad = cem * EM1 + cemp * EMP1 + cd * D1 + cf * F1;
        p += coeff * std::cos(theta);
        pd -= coeff * std::sin(theta) * thetad;
    }
    p = (p + 0.950724f) * DEGREES_2_RADIANS;
    pd *= RATE_CONV_FACTOR;

    // transform parallax to distance (AU, AU/sec)
    const float sp = std::sin(p);
    const float r = EARTH_RADIUS_AU / sp;
    const float rd = -r * pd / sp;

    // longitude, latitude to x,y,z (AU)
    VectorPV<float> v;
    cs2c6({{{el, b}, r}, {{eld, bd}, rd}}, v);

    // mean obliquity
    const float eps = DEGREES_2_RADIANS * (23.45229f - 0.00013f * t);
    const float sin_eps = std::sin(eps);
    const float cos_eps = std::cos(eps);

    // rotate Moon position and velocity into equatorial system
    pv.set_x(v.get_x());
    pv.set_y(v.get_y() * cos_eps - v.get_z() * sin_eps);
    pv.set_z(v.get_y() * sin_eps + v.get_z() * cos_eps);
    pv.set_dx(v.get_dx());
    pv.set_dy(v.get_dy() * cos_eps - v.get_dz() * sin_eps);
    pv.set_dz(v.get_dy() * sin_eps + v.get_dz() * cos_eps);
}

}
