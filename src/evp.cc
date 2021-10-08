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
 * Calculates barycentric and heliocentric velocity and position of the Earth (double precision).
 *
 * This function is accurate enough for many purposes but faster and more compact than the sla::epv() function. The
 * maximum deviations from the JPL DE96 ephemeris are as follows:
 *
 *   barycentric velocity:   0.42 m/s
 *   barycentric position:   6900 km
 *
 *   heliocentric velocity:  0.42 m/s
 *   heliocentric position:  1600 km
 *
 * The function is adapted from the BARVEL and BARCOR subroutines of Stumpff (1980). Most of the changes are merely
 * cosmetic and do not affect the results at all. However, some adjustments have been made so as to give results that
 * refer to the IAU 1976 'FK5' equinox and precession, although the differences these changes make relative to the
 * results from Stumpff's original 'FK4' version are smaller than the inherent accuracy of the algorithm. One minor
 * shortcoming in the original routines that has NOT been corrected is that better numerical accuracy could be
 * achieved if the various polynomial evaluations were nested.
 *
 * Reference:
 *   Stumpff, P., Astron.Astrophys.Suppl.Ser. 41, 1-8 (1980).
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param date TDB (Barycentric Dynamical Time; loosely ET) as a Modified Julian Date (JD-2400000.5).
 * @param deqx Julian Epoch (e.g. 2000.0D0) of mean equator and equinox of the vectors returned; if `deqx` <= 0.0,
 *   then all vectors are referred to the mean equator and equinox (FK5) of epoch `date`.
 * @param bvelo Return value: barycentric velocity (AU/s, Cartesian vector).
 * @param bpos Return value: barycentric position (AU, Cartesian vector).
 * @param hvelo Return value: heliocentric velocity (AU/s, Cartesian vector).
 * @param hpos Return value: heliocentric position (AU, Cartesian vector).
 */
void evp(double date, double deqx,
    Vector<double> bvelo, Vector<double> bpos, Vector<double> hvelo, Vector<double> hpos) {
    float sn[4], forbel[7], sorbel[17], sinlp[4], coslp[4];
    const float& e = sorbel[0];
    const float& g = forbel[0];

    constexpr double DC2PI = 6.2831853071796;
    constexpr float CC2PI = 6.283185f;
    constexpr double DS2R = 0.7272205216643e-4;
    constexpr double B1950 = 1949.9997904423;

    // constants DCFEL[k][i] of fast changing elements
    //   i==0           i==1                  i==2
    static const double DCFEL[8][3] = {
        {1.7400353e+00, 6.2833195099091e+02,  5.2796e-06},
        {6.2565836e+00, 6.2830194572674e+02, -2.6180e-06},
        {4.7199666e+00, 8.3997091449254e+03, -1.9780e-05},
        {1.9636505e-01, 8.4334662911720e+03, -5.6044e-05},
        {4.1547339e+00, 5.2993466764997e+01,  5.8845e-06},
        {4.6524223e+00, 2.1354275911213e+01,  5.6797e-06},
        {4.2620486e+00, 7.5025342197656e+00,  5.5317e-06},
        {1.4740694e+00, 3.8377331909193e+00,  5.6093e-06}
    };
    //
    // constants DCEPS[i] and CCSEL[k][i] of slowly changing elements
    //  i==0            i==1           i==2
    static const double DCEPS[3] = {
        4.093198e-01,  -2.271110e-04, -2.860401e-08
    };
    static const float CCSEL[17][3] = {
        {1.675104e-02f, -4.179579e-05f, -1.260516e-07f},
        {2.220221e-01f,  2.809917e-02f,  1.852532e-05f},
        {1.589963e+00f,  3.418075e-02f,  1.430200e-05f},
        {2.994089e+00f,  2.590824e-02f,  4.155840e-06f},
        {8.155457e-01f,  2.486352e-02f,  6.836840e-06f},
        {1.735614e+00f,  1.763719e-02f,  6.370440e-06f},
        {1.968564e+00f,  1.524020e-02f, -2.517152e-06f},
        {1.282417e+00f,  8.703393e-03f,  2.289292e-05f},
        {2.280820e+00f,  1.918010e-02f,  4.484520e-06f},
        {4.833473e-02f,  1.641773e-04f, -4.654200e-07f},
        {5.589232e-02f, -3.455092e-04f, -7.388560e-07f},
        {4.634443e-02f, -2.658234e-05f,  7.757000e-08f},
        {8.997041e-03f,  6.329728e-06f, -1.939256e-09f},
        {2.284178e-02f, -9.941590e-05f,  6.787400e-08f},
        {4.350267e-02f, -6.839749e-05f, -2.714956e-07f},
        {1.348204e-02f,  1.091504e-05f,  6.903760e-07f},
        {3.106570e-02f, -1.665665e-04f, -1.590188e-07f}
    };
    //
    // constants of the arguments of the short-period perturbations
    // by the planets: DCARGS[k][i]
    //   i==0            i==1
    static const double DCARGS[15][2] = {
        {5.0974222e+00, -7.8604195454652e+02},
        {3.9584962e+00, -5.7533848094674e+02},
        {1.6338070e+00, -1.1506769618935e+03},
        {2.5487111e+00, -3.9302097727326e+02},
        {4.9255514e+00, -5.8849265665348e+02},
        {1.3363463e+00, -5.5076098609303e+02},
        {1.6072053e+00, -5.2237501616674e+02},
        {1.3629480e+00, -1.1790629318198e+03},
        {5.5657014e+00, -1.0977134971135e+03},
        {5.0708205e+00, -1.5774000881978e+02},
        {3.9318944e+00,  5.2963464780000e+01},
        {4.8989497e+00,  3.9809289073258e+01},
        {1.3097446e+00,  7.7540959633708e+01},
        {3.5147141e+00,  7.9618578146517e+01},
        {3.5413158e+00, -5.4868336758022e+02}
    };
    //
    // amplitudes CCAMPS[k][n] of the short-period perturbations
    //  n==0         n==1         n==2         n==3         n==4
    static const float CCAMPS[15][5] = {
        {-2.279594e-5f,  1.407414e-5f,  8.273188e-6f,  1.340565e-5f, -2.490817e-7f},
        {-3.494537e-5f,  2.860401e-7f,  1.289448e-7f,  1.627237e-5f, -1.823138e-7f},
        { 6.593466e-7f,  1.322572e-5f,  9.258695e-6f, -4.674248e-7f, -3.646275e-7f},
        { 1.140767e-5f, -2.049792e-5f, -4.747930e-6f, -2.638763e-6f, -1.245408e-7f},
        { 9.516893e-6f, -2.748894e-6f, -1.319381e-6f, -4.549908e-6f, -1.864821e-7f},
        { 7.310990e-6f, -1.924710e-6f, -8.772849e-7f, -3.334143e-6f, -1.745256e-7f},
        {-2.603449e-6f,  7.359472e-6f,  3.168357e-6f,  1.119056e-6f, -1.655307e-7f},
        {-3.228859e-6f,  1.308997e-7f,  1.013137e-7f,  2.403899e-6f, -3.736225e-7f},
        { 3.442177e-7f,  2.671323e-6f,  1.832858e-6f, -2.394688e-7f, -3.478444e-7f},
        { 8.702406e-6f, -8.421214e-6f, -1.372341e-6f, -1.455234e-6f, -4.998479e-8f},
        {-1.488378e-6f, -1.251789e-5f,  5.226868e-7f, -2.049301e-7f,  0.0f},
        {-8.043059e-6f, -2.991300e-6f,  1.473654e-7f, -3.154542e-7f,  0.0f},
        { 3.699128e-6f, -3.316126e-6f,  2.901257e-7f,  3.407826e-7f,  0.0f},
        { 2.550120e-6f, -1.241123e-6f,  9.901116e-8f,  2.210482e-7f,  0.0f},
        {-6.351059e-7f,  2.341650e-6f,  1.061492e-6f,  2.878231e-7f,  0.0f}
    };
    //
    // constants of the secular perturbations in longitude CCSEC3 and CCSEC[k][n]
    //   n==0          n==1          n==2
    constexpr float CCSEC3 = -7.757020e-08f;
    static const float CCSEC[4][3] = {
        {1.289600e-06f, 5.550147e-01f, 2.076942e+00f},
        {3.102810e-05f, 4.035027e+00f, 3.525565e-01f},
        {9.124190e-06f, 9.990265e-01f, 2.622706e+00f},
        {9.793240e-07f, 5.508259e+00f, 1.559103e+01f}
    };

    // Sidereal rate DCSLD in longitude, rate CCSGD in mean anomaly
    constexpr double DCSLD = 1.990987e-07;
    constexpr float CCSGD = 1.990969e-07f;

    // Some constants used in the calculation of the lunar contribution
    constexpr float CCKM = 3.122140e-05f;
    constexpr float CCMLD = 2.661699e-06f;
    constexpr float CCFDI = 2.399485e-07f;

    // constants DCARGM[k][i] of the arguments of the perturbations of the motion of the Moon
    //   i==0            i==1
    static const double DCARGM[3][2] = {
        {5.1679830e+00,  8.3286911095275e+03},
        {5.4913150e+00, -7.2140632838100e+03},
        {5.9598530e+00,  1.5542754389685e+04}
    };

    //
    // Amplitudes CCAMPM[k][n] of the perturbations of the Moon
    //    n==0          n==1          n==2           n==3
    static const float CCAMPM[3][4] = {
        { 1.097594e-01, 2.896773e-07, 5.450474e-02,  1.438491e-07},
        {-2.223581e-02, 5.083103e-08, 1.002548e-02, -2.291823e-08},
        { 1.148966e-02, 5.658888e-08, 8.249439e-03,  4.063015e-08}
    };

    // CCPAMV[k]=a*M*DL/dt (planets), DC1MME=1-MASS(Earth+Moon)
    static const float CCPAMV[4] = {
        8.326827e-11,1.843484e-11,1.988712e-12,1.881276e-12
    };
    constexpr double DC1MME = 0.99999696;

    // CCPAM[k]=a*M(planets), CCIM=INCLINATION(Moon)
    static const float CCPAM[4] = {
        4.960906e-3, 2.727436e-3, 8.392311e-4, 1.556861e-3
    };
    constexpr float CCIM = 8.978749e-2f;

    //
    // EXECUTION
    // ---------

    // control parameter deq, and time arguments
    bool deq = false;
    if (deqx > 0.0) {
        deq = true;
    }
    const double dt = (date - 15019.5) / 36525.0;
    auto t = (const float) dt;
    const double dtsq = dt * dt;
    auto tsq = (const float) dtsq;

    // values of all elements for the instant date
    int k;
    double dml;
    for (k = 0; k < 8; k++) {
        const double local = std::fmod(DCFEL[k][0] + dt * DCFEL[k][1] + dtsq * DCFEL[k][2], DC2PI);
        if (k == 0) {
            dml = local;
        } else {
            forbel[k - 1] = (float) local;
        }
    }
    const double deps = std::fmod(DCEPS[0] + dt * DCEPS[1] + dtsq * DCEPS[2], DC2PI);
    for (k = 0; k < 17; k++) {
        sorbel[k] = std::fmod(CCSEL[k][0] + t * CCSEL[k][1] + tsq * CCSEL[k][2], CC2PI);
    }

    // secular perturbations in longitude
    float a;
    for (k = 0; k < 4; k++) {
        a = std::fmod(CCSEC[k][1] + t * CCSEC[k][2], CC2PI);
        sn[k] = std::sin(a);
    }

    // periodic perturbations of the EMB (Earth-Moon barycentre)
    float pertl = CCSEC[0][0] * sn[0] + CCSEC[1][0] * sn[1] + (CCSEC[2][0] + t * CCSEC3) * sn[2] + CCSEC[3][0] * sn[3];
    float pertld = 0.0f;
    float pertr = 0.0f;
    float pertrd = 0.0f;
    float cos_a, sin_a;
    for (k = 0; k < 15; k++) {
        a = float(std::fmod(DCARGS[k][0] + dt * DCARGS[k][1], DC2PI));
        cos_a = std::cos(a);
        sin_a = std::sin(a);
        pertl += CCAMPS[k][0] * cos_a + CCAMPS[k][1] * sin_a;
        pertr += CCAMPS[k][2] * cos_a + CCAMPS[k][3] * sin_a;
        if (k < 10) {
            pertld += (CCAMPS[k][1] * cos_a - CCAMPS[k][0] * sin_a) * CCAMPS[k][4];
            pertrd += (CCAMPS[k][3] * cos_a - CCAMPS[k][2] * sin_a) * CCAMPS[k][4];
        }
    }

    // elliptic part of the motion of the EMB
    const float esq = e * e;
    const double dparam = 1.0 - double(esq);
    auto param = (const float) dparam;
    const float two_e = e + e;
    const float two_g = g + g;
    const float phi = two_e * ((1.0f - esq * 0.125f) * std::sin(g) + e * 0.625f * std::sin(two_g) +
        esq * 0.54166667f * std::sin(g + two_g));
    const float f = g + phi;
    const float sin_f = std::sin(f);
    const float cos_f = std::cos(f);
    const double dpsi = dparam / (1.0 + double(e * cos_f));
    const float phid = two_e * CCSGD * ((1.0f + esq * 1.5f) * cos_f + e * (1.25f - sin_f * sin_f * 0.5f));
    const float psid = CCSGD * e * sin_f / std::sqrt(param);

    // perturbed heliocentric motion of the EMB
    const double d1pdro = 1.0 + double(pertr);
    const double drd = d1pdro * (double(psid) + dpsi * double(pertrd));
    const double drld = d1pdro * dpsi * (DCSLD + double(phid) + double(pertld));
    const double dtl = std::fmod(dml + double(phi) + double(pertl), DC2PI);
    const double dsinls = std::sin(dtl);
    const double dcosls = std::cos(dtl);
    double dxhd = drd * dcosls - drld * dsinls;
    double dyhd = drd * dsinls + drld * dcosls;

    // influence of eccentricity, evection and variation on the geocentric motion of the Moon
    pertl = 0.0f;
    pertld = 0.0f;
    float pertp = 0.0f;
    float pertpd = 0.0f;
    for (k = 0; k < 3; k++) {
        a = float(std::fmod(DCARGM[k][0] + dt * DCARGM[k][1], DC2PI));
        sin_a = std::sin(a);
        cos_a = std::cos(a);
        pertl += CCAMPM[k][0] * sin_a;
        pertld += CCAMPM[k][1] * cos_a;
        pertp += CCAMPM[k][2] * cos_a;
        pertpd -= CCAMPM[k][3] * sin_a;
    }

    // heliocentric motion of the Earth
    float tl = forbel[1] + pertl;
    const float sinlm = std::sin(tl);
    const float coslm = std::cos(tl);
    const float sigma = CCKM / (1.0f + pertp);
    a = sigma * (CCMLD + pertld);
    float b = sigma * pertpd;
    dxhd = dxhd + double(a * sinlm) + double(b * coslm);
    dyhd = dyhd - double(a * coslm) + double(b * sinlm);
    const double dzhd = -double(sigma * CCFDI * std::cos(forbel[2]));

    // barycentric motion of the Earth
    double dxbd = dxhd * DC1MME;
    double dybd = dyhd * DC1MME;
    double dzbd = dzhd * DC1MME;
    for (k = 0; k < 4; k++) {
        const float plon = forbel[k + 3];
        const float pomg = sorbel[k + 1];
        const float pecc = sorbel[k + 9];
        tl = std::fmod(plon + 2.0f * pecc * std::sin(plon - pomg), CC2PI);
        sinlp[k] = std::sin(tl);
        coslp[k] = std::cos(tl);
        dxbd += double(CCPAMV[k] * (sinlp[k] + pecc * std::sin(pomg)));
        dybd -= double(CCPAMV[k] * (coslp[k] + pecc * std::cos(pomg)));
        dzbd -= double(CCPAMV[k] * sorbel[k + 13] * std::cos(plon - sorbel[k + 5]));
    }

    // transition to mean equator of date
    const double dcosep = std::cos(deps);
    const double dsinep = std::sin(deps);
    const double dyahd = dcosep * dyhd - dsinep * dzhd;
    const double dzahd = dsinep * dyhd + dcosep * dzhd;
    const double dyabd = dcosep * dybd - dsinep * dzbd;
    const double dzabd = dsinep * dybd + dcosep * dzbd;

    // heliocentric coordinates of the Earth
    const double dr = dpsi * d1pdro;
    const float flatm = CCIM * std::sin(forbel[2]);
    a = sigma * std::cos(flatm);
    const double dxh = dr * dcosls - double(a * coslm);
    const double dyh = dr * dsinls - double(a * sinlm);
    const double dzh = -double(sigma * std::sin(flatm));

    // barycentric coordinates of the Earth
    double dxb = dxh * DC1MME;
    double dyb = dyh * DC1MME;
    double dzb = dzh * DC1MME;
    for (k = 0; k < 4; k++) {
        const float flat = sorbel[k + 13] * std::sin(forbel[k + 3] - sorbel[k + 5]);
        a = CCPAM[k] * (1.0f - sorbel[k + 9] * std::cos(forbel[k + 3] - sorbel[k + 1]));
        b = a * std::cos(flat);
        dxb -= double(b * coslp[k]);
        dyb -= double(b * sinlp[k]);
        dzb -= double(a * std::sin(flat));
    }

    // transition to mean equator of date
    const double dyah = dcosep * dyh - dsinep * dzh;
    const double dzah = dsinep * dyh + dcosep * dzh;
    const double dyab = dcosep * dyb - dsinep * dzb;
    const double dzab = dsinep * dyb + dcosep * dzb;

    // copy result components into vectors, correcting for FK4 equinox
    const double depj = epj(date);
    const double deqcor = DS2R * (0.035 + 0.00085 * (depj - B1950));
    hvelo[0] = dxhd - deqcor * dyahd;
    hvelo[1] = dyahd + deqcor * dxhd;
    hvelo[2] = dzahd;
    bvelo[0] = dxbd - deqcor * dyabd;
    bvelo[1] = dyabd + deqcor * dxbd;
    bvelo[2] = dzabd;
    hpos[0] = dxh - deqcor * dyah;
    hpos[1] = dyah + deqcor * dxh;
    hpos[2] = dzah;
    bpos[0] = dxb - deqcor * dyab;
    bpos[1] = dyab + deqcor * dxb;
    bpos[2] = dzab;

    // was precession to another equinox requested?
    if (deq) {

        // yes: compute precession matrix from MJD date to Julian epoch deqx
        Matrix<double> mat;
        prec(depj, deqx, mat);

        // rotate hvelo
        Vector<double> vec;
        double w;
        int i, j;
        for (j = 0; j < 3; j++) {
            w = 0.0;
            for (i = 0; i < 3; i++) {
                w += mat[j][i] * hvelo[i];
            }
            vec[j] = w;
        }
        for (j = 0; j < 3; j++) {
            hvelo[j] = vec[j];
        }

        // rotate bvelo
        for (j = 0; j < 3; j++) {
            w = 0.0;
            for (i = 0; i < 3; i++) {
                w += mat[j][i] * bvelo[i];
            }
            vec[j] = w;
        }
        for (j = 0; j < 3; j++) {
            bvelo[j] = vec[j];
        }

        // rotate hpos
        for (j = 0; j < 3; j++) {
            w = 0.0;
            for (i = 0; i < 3; i++) {
                w += mat[j][i] * hpos[i];
            }
            vec[j] = w;
        }
        for (j = 0; j < 3; j++) {
            hpos[j] = vec[j];
        }

        // rotate bpos
        for (j = 0; j < 3; j++) {
            w = 0.0;
            for (i = 0; i < 3; i++) {
                w += mat[j][i] * bpos[i];
            }
            vec[j] = w;
        }
        for (j = 0; j < 3; j++) {
            bpos[j] = vec[j];
        }
    }
}

}
