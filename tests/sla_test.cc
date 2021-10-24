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
#include <cstdio>
#include <cstring>
#include <cmath>

#include "sla_test.h"
#include "../src/slalib.h"

namespace sla {

///////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS
///////////////////////////////////////////////////////////////////////////////

// reports error in a function test
static void err(const char* func, const char* test, bool &status) {
    std::printf("Test '%s' of the function %s() FAILED\n", test, func);
    status = false;
}

// validates a character string result
static void vcs(const char* str, const char* str_ok, const char* func, const char* test, bool &status) {
    if (std::strcmp(str, str_ok) != 0) {
        err(func, test, status);
        std::printf("\tExpected: '%s'\n", str_ok);
        std::printf("\t  Actual: '%s'\n", str);
    }
}

// validates an integer result
static void viv(const int val, const int val_ok, const char* func, const char* test, bool &status) {
    if (val != val_ok) {
        err(func, test, status);
        std::printf("\tExpected: %d\n", val_ok);
        std::printf("\t  Actual: %d\n", val);
    }
}

// validates a long result; in FORTRAN implementation of SLALIB, "long" values are INTEGER*4, so C/C++ `int`s
// work just fine, there is no need to use `long`s
static void vlv(const int val, const int val_ok, const char* func, const char* test, bool &status) {
    static_assert(sizeof val >= 4, "`int` values must be at least 32-bit");
    viv(val, val_ok, func, test, status);
}

// validates a double precision floating point result
static void vvd(const double val, const double val_ok, const double tolerance,
    const char* func, const char* test, bool &status) {
    if (std::fabs(val - val_ok) > tolerance)  {
        err(func, test, status);
        std::printf("\tExpected: %30.19f\n", val_ok);
        std::printf("\t  Actual: %30.19f\n", val);
    }
}

///////////////////////////////////////////////////////////////////////////////
// INDIVIDUAL FUNCTION TESTS
///////////////////////////////////////////////////////////////////////////////

// tests sla::airmas() function
static void t_airmas(bool& status) {
    vvd(airmas(1.2354), 3.015698990074724, 1e-12, "sla::airmas", "", status);
}

// tests sla::bear(), sla::dbear(), sla::pav(), and sla::dpav() functions
static void t_bear(bool& status) {
    Vector<float> fv1, fv2;
    Vector<double> dv1, dv2;
    constexpr double a1 = 1.234;
    constexpr double b1 = -0.123;
    constexpr double a2 = 2.345;
    constexpr double b2 = 0.789;

    vvd(bear({float(a1), float(b1)}, {float(a2), float(b2)}),
        0.7045970341781791, 1.0e-6, "sla::bear", "", status);
    vvd(dbear({a1, b1}, {a2, b2}),
        0.7045970341781791, 1.0e-12, "sla::dbear", "", status);
    dcs2c({a1, b1}, dv1);
    dcs2c({a2, b2}, dv2);

    for (int i = 0; i < 3; i++) {
        fv1[i] = float(dv1[i]);
        fv2[i] = float(dv2[i]);
    }
    vvd(pav(fv1, fv2 ), 0.7045970341781791, 1.0e-6, "sla::pav", "", status);
    vvd(dpav(dv1, dv2), 0.7045970341781791, 1.0e-12, "sla::dpav", "", status);
}

// tests sla::caf2r() and sla::daf2r() procedures
static void t_caf2r(bool& status) {
    float sp_radians;
    D2RStatus result = caf2r(76, 54, 32.1f, sp_radians);
    vvd(sp_radians, 1.342313819975276, 1.0e-6, "sla::caf2r", "r", status);
    viv(result, 0, "sla::caf2r", "s", status );

    double dp_radians;
    result = daf2r(76, 54, 32.1, dp_radians);
    vvd(dp_radians, 1.342313819975276, 1.0e-12, "sla::daf2r", "r", status);
    viv(result, 0, "sla::daf2r", "s", status);
}

// tests sla::caldj() procedure
static void t_caldj(bool& status) {
    double mjd;
    auto result = caldj(1999, 12, 31, mjd);
    vvd(mjd, 51543.0, 0.0, "sla::caldj", "d", status);
    viv(result, 0, "sla::caldj", "s", status);
}

// tests sla::clyd() and sla::calyd() procedures
static void t_calyd(bool& status) {
    int j_year, j_day;
    G2JStatus result;

    result = calyd(46, 4, 30, j_year, j_day);
    viv(j_year, 2046, "sla::calyd", "year", status);
    viv(j_day, 120, "sla::calyd", "day", status);
    viv(result, 0, "sla::calyd", "status", status);
    result = clyd (-5000, 1, 1, j_year, j_day);
    viv(result, 1, "sla::clyd", "illegal year", status);
    result = clyd (1900, 0, 1, j_year, j_day);
    viv(result, 2, "sla::clyd", "illegal month", status);
    result = clyd (1900, 2, 29, j_year, j_day);
    viv(j_year, 1900, "sla::clyd", "illegal day (y)", status);
    viv(j_day, 61, "sla::clyd", "illegal day (d)", status);
    viv(result, 3, "sla::clyd", "illegal day (s)", status);
    result = clyd (2000, 2, 29, j_year, j_day);
    viv(j_year, 2000, "sla::clyd", "year", status);
    viv(j_day, 60, "sla::clyd", "day", status);
    viv(result, 0, "sla::clyd", "status", status);
}

// tests sla::clyd() and sla::calyd() procedures
static void t_djcal(bool& status) {
    constexpr double DJM = 50123.9999;
    Date date;
    bool result = djcal( 4, DJM, date);
    viv(date.d_year, 1996, "sla::djcal", "year", status);
    viv(date.d_month, 2, "sla::djcal", "month", status);
    viv(date.d_day, 10, "sla::djcal", "day", status);
    viv(date.d_ifraction, 9999, "sla::djcal", "fraction", status);
    viv((int) result, 0, "sla::djcal", "status", status);

    result = djcl(DJM, date);
    viv(date.d_year, 1996, "sla::djcl", "year", status);
    viv(date.d_month, 2, "sla::djcl", "month", status);
    viv(date.d_day, 10, "sla::djcl", "day", status);
    vvd(date.d_fraction, 0.9999, 1.0e-7, "sla::djcl", "fraction", status);
    viv((int) result, 0, "sla::djcl", "status", status);
}

// tests sla::cc2s() and dcc2s() procedures
static void t_cc2s(bool& status) {
    const Vector<float> v = {100.0f, -50.0f, 25.0f};
    Spherical<float> sp_spherical;
    cc2s(v, sp_spherical);
    vvd(sp_spherical.get_ra(), -0.4636476090008061, 1.0e-6, "sla::cc2s", "ra", status );
    vvd(sp_spherical.get_dec(), 0.2199879773954594, 1.0e-6, "sla::cc2s", "dec", status );

    const Vector<double> dv = {100.0, -50.0, 25.0};
    Spherical<double> dp_spherical;
    dcc2s(dv, dp_spherical);
    vvd(dp_spherical.get_ra(), -0.4636476090008061, 1.0e-12, "sla::dcc2s", "ra", status );
    vvd(dp_spherical.get_dec(), 0.2199879773954594, 1.0e-12, "sla::dcc2s", "dec", status );
}

// tests sla::cldj() procedure
static void t_cldj(bool& status) {
    double mjd;
    auto result = cldj(1899, 12, 31, mjd);
    vvd(mjd, 15019.0, 0.0, "sla::cldj", "d", status);
    viv(result, 0, "sla::cldj", "s", status);
}

// tests sla::e2h(), sla::de2h(), sla::h2e(), and sla::dh2e() procedures
static void t_e2h(bool& status) {
    Spherical<double> d_dir = {-0.3, -1.1};
    const double d_phi = -0.7;
    Spherical<float> f_dir = {(float) d_dir.get_ha(), (float) d_dir.get_dec()};
    const float f_phi = (float) d_phi;

    double d_azimuth, d_elevation;
    de2h(d_dir, d_phi, d_azimuth, d_elevation);
    vvd(d_azimuth, 2.820087515852369, 1.0e-12, "sla::de2h", "azimuth", status);
    vvd(d_elevation, 1.132711866443304, 1.0e-12, "sla::de2h", "elevation", status);

    float azimuth, elevation;
    e2h(f_dir, f_phi, azimuth, elevation);
    vvd(azimuth, 2.820087515852369, 1.0e-6, "sla::e2h", "azimuth", status);
    vvd(elevation, 1.132711866443304, 1.0e-6, "sla::e2h", "elevation", status);

    dh2e(d_azimuth, d_elevation, d_phi, d_dir);
    vvd(d_dir.get_ha(), -0.3, 1.0e-12, "sla::dh2e", "ha", status);
    vvd(d_dir.get_dec(), -1.1, 1.0e-12, "sla::dh2e", "dec", status);

    h2e(azimuth, elevation, f_phi, f_dir);
    vvd(f_dir.get_ha(), -0.3, 1.0e-6, "sla::h2e", "ha", status);
    vvd(f_dir.get_dec(), -1.1, 1.0e-6, "sla::h2e", "dec", status);
}

/*
 * tests all the 3-component vector and 3x3 matrix procedures:
 *
 *   sla::av2m()   sla::dav2m()
 *   sla::cs2c()   sla::dcs2c()      (these two take structures in C++ implementation)
 *   sla::euler()  sla::deuler()
 *   sla::imxv()   sla::dimxv()
 *   sla::m2av()   sla::dm2av()
 *   sla::mxm()    sla::dmxm()
 *   sla::mxv()    sla::dmxv()
 *   sla::vdv()    sla::dvdv()
 *   sla::vn()     sla::dvn()        (these two return values in the C++ implementation)
 *   sla::vxv()    sla::dvxv()
 *
 * sla::cc2s() and sla::dcc2s() are tested in separate procedure t_cc2s(), even though FORTRAN implementation lists
 * them as tested in the T_VECMAT subroutine.
 */
static void t_vecmat(bool& status) {
    // tolerances for [most] single and double precision functions' tests, respectively
    constexpr double sp_tolerance = 1.0e-6;
    constexpr double dp_tolerance = 1.0e-12;

    // make a rotation matrix
    Vector<float> av;
    Matrix<float> rm1;
    av[0] = -0.123f;
    av[1] = 0.0987f;
    av[2] = 0.0654f;
    av2m(av, rm1);
    vvd(rm1[0][0], 0.9930075842721269, sp_tolerance, "sla::av2m", "00", status);
    vvd(rm1[0][1], 0.05902743090199868, sp_tolerance, "sla::av2m", "01", status);
    vvd(rm1[0][2], -0.1022335560329612, sp_tolerance, "sla::av2m", "02", status);
    vvd(rm1[1][0], -0.07113807138648245, sp_tolerance, "sla::av2m", "10", status);
    vvd(rm1[1][1], 0.9903204657727545, sp_tolerance, "sla::av2m", "11", status);
    vvd(rm1[1][2], -0.1191836812279541, sp_tolerance, "sla::av2m", "12", status);
    vvd(rm1[2][0], 0.09420887631983825, sp_tolerance, "sla::av2m", "20", status);
    vvd(rm1[2][1], 0.1256229973879967, sp_tolerance, "sla::av2m", "21", status);
    vvd(rm1[2][2], 0.9875948309655174, sp_tolerance, "sla::av2m", "22", status);

    // make another
    Matrix<float> rm2;
    euler("YZY", 2.345E0, -0.333E0, 2.222E0, rm2);
    vvd(rm2[0][0], -0.1681574770810878, sp_tolerance, "sla::euler", "00", status);
    vvd(rm2[0][1], 0.1981362273264315, sp_tolerance, "sla::euler", "01", status);
    vvd(rm2[0][2], 0.9656423242187410, sp_tolerance, "sla::euler", "02", status);
    vvd(rm2[1][0], -0.2285369373983370, sp_tolerance, "sla::euler", "10", status);
    vvd(rm2[1][1], 0.9450659587140423, sp_tolerance, "sla::euler", "11", status);
    vvd(rm2[1][2], -0.2337117924378156, sp_tolerance, "sla::euler", "12", status);
    vvd(rm2[2][0], -0.9589024617479674, sp_tolerance, "sla::euler", "20", status);
    vvd(rm2[2][1], -0.2599853247796050, sp_tolerance, "sla::euler", "21", status);
    vvd(rm2[2][2], -0.1136384607117296, sp_tolerance, "sla::euler", "22", status);

    // combine them
    Matrix<float> rm;
    mxm(rm2, rm1, rm);
    vvd(rm[0][0], -0.09010460088585805, sp_tolerance, "sla::mxm", "00", status);
    vvd(rm[0][1], 0.3075993402463796, sp_tolerance, "sla::mxm", "01", status);
    vvd(rm[0][2], 0.9472400998581048, sp_tolerance, "sla::mxm", "02", status);
    vvd(rm[1][0], -0.3161868071070688, sp_tolerance, "sla::mxm", "10", status);
    vvd(rm[1][1], 0.8930686362478707, sp_tolerance, "sla::mxm", "11", status);
    vvd(rm[1][2], -0.3200848543149236, sp_tolerance, "sla::mxm", "12", status);
    vvd(rm[2][0], -0.9444083141897035, sp_tolerance, "sla::mxm", "20", status);
    vvd(rm[2][1], -0.3283459407855694, sp_tolerance, "sla::mxm", "21", status);
    vvd(rm[2][2], 0.01678926022795169, sp_tolerance, "sla::mxm", "22", status);

    // create a vector
    Vector<float> v1;
    cs2c({3.0123f, -0.999f}, v1 );
    vvd(v1[0], -0.5366267667260525, sp_tolerance, "sla::cs2c", "x", status);
    vvd(v1[1], 0.06977111097651444, sp_tolerance, "sla::cs2c", "y", status);
    vvd(v1[2], -0.8409302618566215, sp_tolerance, "sla::cs2c", "z", status);

    // rotate the vector using the two matrices sequentially
    Vector<float> v2, v3;
    mxv(rm1, v1, v2);
    mxv(rm2, v2, v3);
    vvd(v3[0], -0.7267487768696160, sp_tolerance, "sla::mxv", "x", status);
    vvd(v3[1], 0.5011537352639822, sp_tolerance, "sla::mxv", "y", status);
    vvd(v3[2], 0.4697671220397141, sp_tolerance, "sla::mxv", "z", status);

    // de-rotate the vector using the combined matrix
    Vector<float> v4;
    imxv(rm, v3, v4);
    vvd(v4[0], -0.5366267667260526, sp_tolerance, "sla::imxv", "x", status);
    vvd(v4[1], 0.06977111097651445, sp_tolerance, "sla::imxv", "y", status);
    vvd(v4[2], -0.8409302618566215, sp_tolerance, "sla::imxv", "z", status);

    // convert the combined matrix into an axial vector
    Vector<float> v5;
    m2av(rm, v5);
    vvd(v5[0], 0.006889040510209034, sp_tolerance, "sla::m2av", "x", status);
    vvd(v5[1], -1.577473205461961, sp_tolerance, "sla::m2av", "y", status);
    vvd(v5[2], 0.5201843672856759, sp_tolerance, "sla::m2av", "z", status);

    // multiply the axial vector by a scalar and then normalize
    for (int i = 0; i < 3; i++) {
        v5[i] *= 1000.0f;
    }
    Vector<float> v6;
    float vm = vn(v5, v6);
    vvd(v6[0], 0.004147420704640065, sp_tolerance, "sla::vn", "x", status);
    vvd(v6[1], -0.9496888606842218, sp_tolerance, "sla::vn", "y", status);
    vvd(v6[2], 0.3131674740355448, sp_tolerance, "sla::vn", "z", status);
    vvd(vm, 1661.042127339937, 1.0e-3, "sla::vn", "m", status);

    // calculate dot product with the original vector
    vvd(vdv(v6, v1), -0.3318384698006295, sp_tolerance, "sla::vn", "", status);

    // calculate cross product with the original vector
    Vector<float> v7;
    vxv(v6, v1, v7);
    vvd(v7[0], 0.7767720597123304, sp_tolerance, "sla::vxv", "x", status);
    vvd(v7[1], -0.1645663574562769, sp_tolerance, "sla::vxv", "y", status);
    vvd(v7[2], -0.5093390925544726, sp_tolerance, "sla::vxv", "z", status);

    // do same tests in double precision
    Vector<double> dav;
    Matrix<double> drm1;
    dav[0] = -0.123;
    dav[1] = 0.0987;
    dav[2] = 0.0654;
    dav2m(dav, drm1);
    vvd(drm1[0][0], 0.9930075842721269, dp_tolerance, "sla::dav2m", "00", status);
    vvd(drm1[0][1], 0.05902743090199868, dp_tolerance, "sla::dav2m", "01", status);
    vvd(drm1[0][2], -0.1022335560329612, dp_tolerance, "sla::dav2m", "02", status);
    vvd(drm1[1][0], -0.07113807138648245, dp_tolerance, "sla::dav2m", "10", status);
    vvd(drm1[1][1], 0.9903204657727545, dp_tolerance, "sla::dav2m", "11", status);
    vvd(drm1[1][2], -0.1191836812279541, dp_tolerance, "sla::dav2m", "12", status);
    vvd(drm1[2][0], 0.09420887631983825, dp_tolerance, "sla::dav2m", "20", status);
    vvd(drm1[2][1], 0.1256229973879967, dp_tolerance, "sla::dav2m", "21", status);
    vvd(drm1[2][2], 0.9875948309655174, dp_tolerance, "sla::dav2m", "22", status);

    Matrix<double> drm2;
    deuler("YZY", 2.345, -0.333, 2.222, drm2);
    vvd(drm2[0][0], -0.1681574770810878, dp_tolerance, "sla::deuler", "00", status);
    vvd(drm2[0][1], 0.1981362273264315, dp_tolerance, "sla::deuler", "01", status);
    vvd(drm2[0][2], 0.9656423242187410, dp_tolerance, "sla::deuler", "02", status);
    vvd(drm2[1][0], -0.2285369373983370, dp_tolerance, "sla::deuler", "10", status);
    vvd(drm2[1][1], 0.9450659587140423, dp_tolerance, "sla::deuler", "11", status);
    vvd(drm2[1][2], -0.2337117924378156, dp_tolerance, "sla::deuler", "12", status);
    vvd(drm2[2][0], -0.9589024617479674, dp_tolerance, "sla::deuler", "20", status);
    vvd(drm2[2][1], -0.2599853247796050, dp_tolerance, "sla::deuler", "21", status);
    vvd(drm2[2][2], -0.1136384607117296, dp_tolerance, "sla::deuler", "22", status);

    Matrix<double> drm;
    dmxm(drm2, drm1, drm);
    vvd(drm[0][0], -0.09010460088585805, dp_tolerance, "sla::dmxm", "00", status);
    vvd(drm[0][1], 0.3075993402463796, dp_tolerance, "sla::dmxm", "01", status);
    vvd(drm[0][2], 0.9472400998581048, dp_tolerance, "sla::dmxm", "02", status);
    vvd(drm[1][0], -0.3161868071070688, dp_tolerance, "sla::dmxm", "10", status);
    vvd(drm[1][1], 0.8930686362478707, dp_tolerance, "sla::dmxm", "11", status);
    vvd(drm[1][2], -0.3200848543149236, dp_tolerance, "sla::dmxm", "12", status);
    vvd(drm[2][0], -0.9444083141897035, dp_tolerance, "sla::dmxm", "20", status);
    vvd(drm[2][1], -0.3283459407855694, dp_tolerance, "sla::dmxm", "21", status);
    vvd(drm[2][2], 0.01678926022795169, dp_tolerance, "sla::dmxm", "22", status);

    Vector<double> dv1;
    dcs2c({3.0123, -0.999}, dv1);
    vvd(dv1[0], -0.5366267667260525, dp_tolerance, "sla::dcs2c", "x", status);
    vvd(dv1[1], 0.06977111097651444, dp_tolerance, "sla::dcs2c", "y", status);
    vvd(dv1[2], -0.8409302618566215, dp_tolerance, "sla::dcs2c", "z", status);

    Vector<double> dv2, dv3;
    dmxv(drm1, dv1, dv2);
    dmxv(drm2, dv2, dv3);
    vvd(dv3[0], -0.7267487768696160, dp_tolerance, "sla::dmxv", "x", status);
    vvd(dv3[1], 0.5011537352639822, dp_tolerance, "sla::dmxv", "y", status);
    vvd(dv3[2], 0.4697671220397141, dp_tolerance, "sla::dmxv", "z", status);

    Vector<double> dv4;
    dimxv(drm, dv3, dv4);
    vvd(dv4[0], -0.5366267667260526, dp_tolerance, "sla::dimxv", "x", status);
    vvd(dv4[1], 0.06977111097651445, dp_tolerance, "sla::dimxv", "y", status);
    vvd(dv4[2], -0.8409302618566215, dp_tolerance, "sla::dimxv", "z", status);

    Vector<double> dv5;
    dm2av(drm, dv5);
    vvd(dv5[0], 0.006889040510209034, dp_tolerance, "sla::dm2av", "x", status);
    vvd(dv5[1], -1.577473205461961, dp_tolerance, "sla::dm2av", "y", status);
    vvd(dv5[2], 0.5201843672856759, dp_tolerance, "sla::dm2av", "z", status);

    for (int j = 0; j < 3; j++) {
        dv5[j] *= 1000.0;
    }
    Vector<double> dv6;
    double dvm = dvn(dv5, dv6);
    vvd(dv6[0], 0.004147420704640065, dp_tolerance, "sla::dvn", "x", status);
    vvd(dv6[1], -0.9496888606842218, dp_tolerance, "sla::dvn", "y", status);
    vvd(dv6[2], 0.3131674740355448, dp_tolerance, "sla::dvn", "z", status);
    vvd(dvm, 1661.042127339937, 1.0e-9, "sla::dvn", "M", status);

    vvd(dvdv(dv6, dv1), -0.3318384698006295, dp_tolerance, "sla::dvn", "", status);

    Vector<double> dv7;
    dvxv(dv6, dv1, dv7);
    vvd(dv7[0], 0.7767720597123304, dp_tolerance, "sla::dvxv", "x", status);
    vvd(dv7[1], -0.1645663574562769, dp_tolerance, "sla::dvxv", "y", status);
    vvd(dv7[2], -0.5093390925544726, dp_tolerance, "sla::dvxv", "z", status);
}

// tests sla::zd() function
static void t_zd(bool& status) {
    vvd(zd({-1.023, -0.876}, -0.432), 0.8963914139430839, 1.0e-12, "sla::zd", "", status);
}

// tests sla::cd2tf() and sla::dd2tf() procedures
static void t_cd2tf(bool& status) {
    ConversionResult result;
    cd2tf(4, -0.987654321, result);
    viv((int) result.get_sign(), int('-'), "sla::cd2tf", "sign", status);
    viv(result.get_hours(), 23, "sla::cd2tf", "hours", status);
    viv(result.get_minutes(), 42, "sla::cd2tf", "minutes", status);
    viv(result.get_seconds(), 13, "sla::cd2tf", "seconds", status);
    vvd(double(result.get_fraction()), 3333.0, 1000.0, "sla::cd2tf", "fraction", status);

    dd2tf(4, -0.987654321, result);
    viv((int) result.get_sign(), int('-'), "sla::dd2tf", "sign", status);
    viv(result.get_hours(), 23, "sla::dd2tf", "hours", status);
    viv(result.get_minutes(), 42, "sla::dd2tf", "minutes", status);
    viv(result.get_seconds(), 13, "sla::dd2tf", "seconds", status);
    viv(result.get_fraction(), 3333, "sla::dd2tf", "fraction", status);
}

// tests sla::cr2af() and sla::dr2af() procedures
static void t_cr2af(bool& status) {
    ConversionResult result;
    cr2af(4, 2.345f, result);
    viv((int) result.get_sign(), int('+'), "sla::cr2af", "sign", status);
    viv(result.get_degrees(), 134, "sla::cr2af", "degrees", status);
    viv(result.get_arcminutes(), 21, "sla::cr2af", "arcminutes", status);
    viv(result.get_arcseconds(), 30, "sla::cr2af", "arcseconds", status);
    vvd(double(result.get_fraction()), 9706.0, 1000.0, "sla::cr2af", "fraction", status);

    dr2af(4, 2.345, result);
    viv((int) result.get_sign(), int('+'), "sla::dr2af", "sign", status);
    viv(result.get_hours(), 134, "sla::dr2af", "degrees", status);
    viv(result.get_minutes(), 21, "sla::dr2af", "arcminutes", status);
    viv(result.get_seconds(), 30, "sla::dr2af", "arcseconds", status);
    viv(result.get_fraction(), 9706, "sla::dr2af", "fraction", status);
}

// tests sla::cr2tf() and sla::dr2tf() procedures
static void t_cr2tf(bool& status) {
    ConversionResult result;
    cr2tf(4, -3.01234, result);
    viv((int) result.get_sign(), int('-'), "sla::cr2tf", "sign", status);
    viv(result.get_hours(), 11, "sla::cr2tf", "hours", status);
    viv(result.get_minutes(), 30, "sla::cr2tf", "minutes", status);
    viv(result.get_seconds(), 22, "sla::cr2tf", "seconds", status);
    vvd(double(result.get_fraction()), 6484.0, 1000.0, "sla::cr2tf", "fraction", status);

    dr2tf( 4, -3.01234, result);
    viv((int) result.get_sign(), int('-'), "sla::dr2tf", "sign", status);
    viv(result.get_hours(), 11, "sla::dr2tf", "hours", status);
    viv(result.get_minutes(), 30, "sla::dr2tf", "minutes", status);
    viv(result.get_seconds(), 22, "sla::dr2tf", "seconds", status);
    viv(result.get_fraction(), 6484, "sla::dr2tf", "fraction", status);
}

// tests sla::ctf2d() and sla::dtf2d() procedures
static void t_ctf2d(bool& status) {
    float sp_days;
    T2DStatus result = ctf2d (23, 56, 59.1, sp_days);
    vvd(double(sp_days), 0.99790625, 1.0e-6, "sla::ctf2d", "days", status);
    viv(result, T2D_OK, "sla::ctf2d", "result", status);

    double dp_days;
    result = dtf2d (23, 56, 59.1, dp_days);
    vvd(dp_days, 0.99790625, 1.0e-12, "sla::dtf2d", "days", status);
    viv(result, T2D_OK, "sla::dtf2d", "result", status);
}

// tests sla::ctf2r() and sla::dtf2r() procedures
static void t_ctf2r(bool& status) {
    float sp_radians;
    T2DStatus result = ctf2r(23, 56, 59.1, sp_radians);
    vvd (double(sp_radians), 6.270029887942679, 1.0e-6, "sla::ctf2r", "radians", status);
    viv (result, T2D_OK, "sla::ctf2r", "result", status);

    double dp_radians;
    result = dtf2r(23, 56, 59.1, dp_radians);
    vvd (dp_radians, 6.270029887942679, 1.0e-12, "sla::dtf2r", "radians", status);
    viv (result, T2D_OK, "sla::dtf2r", "result", status);
}

// tests sla::altaz() function
static void t_dat(bool& status) {
    vvd(dat(43900.0), 18.0, 0.0, "sla::dat", "43900", status);
    vvd(dtt(40404.0), 39.709746, 1.0e-12, "sla::dtt", "40404", status);
    vvd(dt(500.0), 4686.7, 1.0e-10, "sla::dt", "500", status);
    vvd(dt(1400.0), 408.0, 1.0e-11, "sla::dt", "1400", status);
    vvd(dt(1950.0), 27.99145626, 1.0e-12, "sla::dt", "1950", status);
}


// tests sla::range() and sla::drange() functions
static void t_range(bool& status) {
    vvd(range(-4.0f), 2.283185307179586, 1.0e-6, "sla::range", "float", status);
    vvd(drange(-4.0), 2.283185307179586, 1.0e-12, "sla::drange", "double", status);
}

// tests sla::ranorm() and sla::dranrm() functions
static void t_ranorm(bool& status) {
    vvd(ranorm(-0.1), 6.183185307179587, 1.0e-5, "sla::ranorm", "float", status);
    vvd(dranrm(-0.1), 6.183185307179587, 1.0e-12, "sla::dranrm", "double", status);
}

// tests sla::refro(), sla::refcoq(), sla::refco(), sla::atmdsp(), sla::dcs2c(), sla::refv(), and sla::refz() functions
static void t_ref(bool& status) {
    double ref = refro(1.4, 3456.7, 280.0, 678.9, 0.9, 0.55, -0.3, 0.006, 1.0e-9);
    vvd(ref, 0.00106715763018568, 1.0e-12, "sla::refro", "optical", status);

    ref = refro(1.4, 3456.7, 280.0, 678.9, 0.9, 1000.0, -0.3, 0.006, 1.0e-9);
    vvd(ref, 0.001296416185295403, 1.0e-12, "sla::refro", "radio", status);

    double refa, refb;
    refcoq(275.9, 709.3, 0.9, 101.0, refa, refb);
    vvd(refa, 2.324736903790639e-4, 1.0e-12, "sla::refcoq", "a/r", status);
    vvd(refb, -2.442884551059e-7, 1.0e-15, "sla::refcoq", "b/r", status);

    refco(2111.1, 275.9, 709.3, 0.9, 101.0, -1.03, 0.0067, 1.0e-12, refa, refb);
    vvd(refa, 2.324673985217244e-4, 1.0e-12, "sla::refco", "a/r", status);
    vvd(refb, -2.265040682496e-7, 1.0e-15, "sla::refco", "b/r", status);

    refcoq(275.9, 709.3, 0.9, 0.77, refa, refb);
    vvd(refa, 2.007406521596588e-4, 1.0e-12, "sla::refcoq", "a", status);
    vvd(refb, -2.264210092590e-7, 1.0e-15, "sla::refcoq", "b", status);

    refco(2111.1, 275.9, 709.3, 0.9, 0.77, -1.03, 0.0067, 1.0e-12, refa, refb);
    vvd(refa, 2.007202720084551e-4, 1.0e-12, "sla::refco", "a", status);
    vvd(refb, -2.223037748876e-7, 1.0e-15, "sla::refco", "b", status);

    double refa2, refb2;
    atmdsp(275.9, 709.3, 0.9, 0.77, refa, refb, 0.5, refa2, refb2);
    vvd(refa2, 2.034523658888048e-4, 1.0e-12, "sla::atmdsp", "a", status);
    vvd(refb2, -2.250855362179e-7, 1.0e-15, "sla::atmdsp", "b", status);

    Vector<double> cartesian1, cartesian2;
    dcs2c({0.345, 0.456}, cartesian1);
    refv(cartesian1, refa, refb, cartesian2);
    vvd(cartesian2[0], 0.8447487047790478, 1.0e-12, "sla::refv", "x1", status);
    vvd(cartesian2[1], 0.3035794890562339, 1.0e-12, "sla::refv", "y1", status);
    vvd(cartesian2[2], 0.4407256738589851, 1.0e-12, "sla::refv", "z1", status);

    dcs2c({3.7, 0.03}, cartesian1);
    refv(cartesian1, refa, refb, cartesian2);
    vvd(cartesian2[0], -0.8476187691681673, 1.0e-12, "sla::refv", "x2", status);
    vvd(cartesian2[1], -0.5295354802804889, 1.0e-12, "sla::refv", "y2", status);
    vvd(cartesian2[2], 0.0322914582168426, 1.0e-12, "sla::refv", "z2", status);

    double zr = refz(0.567, refa, refb);
    vvd(zr, 0.566872285910534, 1.0e-12, "sla::refz", "hi el", status);

    zr = refz(1.55, refa, refb);
    vvd(zr, 1.545697350690958, 1.0e-12, "sla::refz", "lo el", status);
}

// tests sla::ecmat() function
static void t_ecmat(bool& status) {
    Matrix<double> rm;
    ecmat(41234.0, rm);
    vvd(rm[0][0], 1.0, 1.0e-12, "sla::ecmat", "00", status);
    vvd(rm[0][1], 0.0, 1.0e-12, "sla::ecmat", "01", status);
    vvd(rm[0][2], 0.0, 1.0e-12, "sla::ecmat", "02", status);
    vvd(rm[1][0], 0.0, 1.0e-12, "sla::ecmat", "10", status);
    vvd(rm[1][1], 0.917456575085716, 1.0e-12, "sla::ecmat", "11", status);
    vvd(rm[1][2], 0.397835937079581, 1.0e-12, "sla::ecmat", "12", status);
    vvd(rm[2][0], 0.0, 1.0e-12, "sla::ecmat", "20", status);
    vvd(rm[2][1], -0.397835937079581, 1.0e-12, "sla::ecmat", "21", status);
    vvd(rm[2][2], 0.917456575085716, 1.0e-12, "sla::ecmat", "22", status);
}

// tests sla::dmat() function
static void t_dmat(bool& status) {
    Matrix<double> mat = {
        {2.22,     1.6578,     1.380522    },
        {1.6578,   1.380522,   1.22548578  },
        {1.380522, 1.22548578, 1.1356276122}
    };
    Vector<double> vec = {2.28625, 1.7128825, 1.429432225};
    double det;
    int ws[3];
    bool singular = dmat(mat, vec, det, ws);

    vvd(mat[0][0], 18.02550629769198, 1.0e-10, "std::dmat", "00", status);
    vvd(mat[0][1], -52.16386644917280607, 1.0e-10, "std::dmat", "01", status);
    vvd(mat[0][2], 34.37875949717850495, 1.0e-10, "std::dmat", "02", status);
    vvd(mat[1][0], -52.16386644917280607, 1.0e-10, "std::dmat", "10", status);
    vvd(mat[1][1], 168.1778099099805627, 1.0e-10, "std::dmat", "11", status);
    vvd(mat[1][2], -118.0722869694232670, 1.0e-10, "std::dmat", "12", status);
    vvd(mat[2][0], 34.37875949717850495, 1.0e-10, "std::dmat", "20", status);
    vvd(mat[2][1], -118.0722869694232670, 1.0e-10, "std::dmat", "21", status);
    vvd(mat[2][2], 86.50307003740151262, 1.0e-10, "std::dmat", "22", status);
    vvd(vec[0], 1.002346480763383, 1.0e-12, "std::dmat", "v0", status);
    vvd(vec[1], 0.03285594016974583489, 1.0e-12, "std::dmat", "v1", status);
    vvd(vec[2], 0.004760688414885247309, 1.0e-12, "std::dmat", "v2", status);
    vvd(det, 0.003658344147359863, 1.0e-12, "std::dmat", "d", status);
    viv((int) singular, 0, "std::dmat", "singular", status);
}

// tests sla::smat() function
static void t_smat(bool& status) {
    Matrix<float> a = {
        {2.22f, 1.6578f, 1.380522f},
        {1.6578f, 1.380522f, 1.22548578f},
        {1.380522f, 1.22548578f, 1.1356276122f}
    };
    Vector<float> v = {2.28625,  1.7128825,  1.429432225};
    float d;
    int w[3];
    bool singular = smat(3, (float*) a, v, d, w);

    vvd(a[0][0], 18.02550629769198, 1.0e-2, "sla::smat", "00", status);
    vvd(a[0][1], -52.16386644917481, 1.0e-2, "sla::smat", "01", status);
    vvd(a[0][2], 34.37875949717994, 1.0e-2, "sla::smat", "02", status);
    vvd(a[1][0], -52.16386644917477, 1.0e-2, "sla::smat", "10", status);
    vvd(a[1][1], 168.1778099099869, 1.0e-1, "sla::smat", "11", status);
    vvd(a[1][2], -118.0722869694278, 1.0e-2, "sla::smat", "12", status);
    vvd(a[2][0], 34.37875949717988, 1.0e-2, "sla::smat", "20", status);
    vvd(a[2][1], -118.07228696942770, 1.0e-2, "sla::smat", "21", status);
    vvd(a[2][2], 86.50307003740468, 1.0e-2, "sla::smat", "22", status);
    vvd(v[0], 1.002346480763383, 1.0e-4, "sla::smat", "v0", status);
    vvd(v[1], 0.0328559401697292, 1.0e-4, "sla::smat", "v1", status);
    vvd(v[2], 0.004760688414898454, 1.0e-4, "sla::smat", "v2", status);
    vvd(d, 0.003658344147359863, 1.0e-4, "sla::smat", "d", status);
    viv((int) singular, 0, "sla::smat", "singular", status);
}

// tests sla::altaz() function
static void t_altaz(bool& status) {
    AltazMount am;
    altaz({0.7, -0.7}, -0.65, am);

    vvd(am.get_azimuth(), 4.400560746660174, 1.0e-12, "sla::altaz", "azimuth", status);
    vvd(am.get_az_velocity(), -0.2015438937145421, 1.0e-13, "sla::altaz", "az_vel", status);
    vvd(am.get_az_acceleration(), -0.4381266949668748, 1.0e-13, "sla::altaz", "az_accel", status);
    vvd(am.get_elevation(), 1.026646506651396, 1.0e-12, "sla::altaz", "elevation", status);
    vvd(am.get_el_velocity(), -0.7576920683826450, 1.0e-13, "sla::altaz", "el_vel", status);
    vvd(am.get_el_acceleration(), 0.04922465406857453, 1.0e-14, "sla::altaz", "el_accel", status);
    vvd(am.get_pangle(), 1.707639969653937, 1.0e-12, "sla::altaz", "pangle", status);
    vvd(am.get_pa_velocity(), 0.4717832355365627, 1.0e-13, "sla::altaz", "pa_vel", status);
    vvd(am.get_pa_acceleration(), -0.2957914128185515, 1.0e-13, "sla::altaz", "pa_accel", status);
}

// tests sla::nut(), sla::nutc(), and sla::nutc80() functions
static void t_nut(bool& status) {
    Matrix<double> mat;
    nut(46012.34, mat);
    vvd(mat[0][0],  9.999999969492166e-1, 1.0e-12, "sla::nut", "00", status);
    vvd(mat[0][1],  7.166577986249302e-5, 1.0e-12, "sla::nut", "01", status);
    vvd(mat[0][2],  3.107382973077677e-5, 1.0e-12, "sla::nut", "02", status);
    vvd(mat[1][0], -7.166503970900504e-5, 1.0e-12, "sla::nut", "10", status);
    vvd(mat[1][1],  9.999999971483732e-1, 1.0e-12, "sla::nut", "11", status);
    vvd(mat[1][2], -2.381965032461830e-5, 1.0e-12, "sla::nut", "12", status);
    vvd(mat[2][0], -3.107553669598237e-5, 1.0e-12, "sla::nut", "20", status);
    vvd(mat[2][1],  2.381742334472628e-5, 1.0e-12, "sla::nut", "21", status);
    vvd(mat[2][2],  9.999999992335206818e-1, 1.0e-12, "sla::nut", "22", status);

    double psi, eps, eps0;
    nutc(50123.4, psi, eps, eps0);
    vvd(psi, 3.523550954747999709e-5, 1.0e-17, "sla::nutc", "psi", status);
    vvd(eps, -4.143371566683342e-5, 1.0e-17, "sla::nutc", "eps", status);
    vvd(eps0, 0.4091014592901651, 1.0e-12, "sla::nutc", "eps0", status);

    nutc80(50123.4, psi, eps, eps0);
    vvd(psi, 3.537714281665945321e-5, 1.0e-17, "sla::nutc80", "psi", status);
    vvd(eps, -4.140590085987148317e-5, 1.0e-17, "sla::nutc80", "deps", status);
    vvd(eps0, 0.4091016349007751, 1.0e-12, "sla::nutc80", "eps0", status);
}

// tests sla::epj2d() function
static void t_epj2d(bool& status) {
    vvd(epj2d(2010.077), 55225.124250, 1.0e-6, "sla::epj2d", "", status);
}

// tests sla::epj() function
static void t_epj(bool& status) {
    vvd(epj(42999.0 ), 1976.603696098563, 1.0e-7, "sla::epj", "", status);
}

// tests sla::epb2d() function
static void t_epb2d(bool& status) {
    vvd(epb2d(1975.5), 42595.5995279655, 1.0e-7, "sla::epb2d", "", status);
}

// tests sla::epb() function
static void t_epb(bool& status) {
    vvd(epb(45123.0), 1982.419793168669, 1.0e-8, "sla::epb", "", status);
}

// tests sla::epco() function
static void t_epco(bool& status) {
    vvd(epco('B', 'J', 2000.0 ), 2000.001277513665, 1.0e-7, "sla::epco", "bj", status);
    vvd(epco('J', 'B', 1950.0 ), 1949.999790442300, 1.0e-7, "sla::epco", "jb", status);
    vvd(epco('J', 'J', 2000.0 ), 2000, 1.0e-7, "sla::epco", "jj", status);
}

// tests sla::prec() and sla::precl() functions
static void t_prec(bool& status) {
    Matrix<double> mat;

    prec(1925.0, 1975.0, mat);
    vvd(mat[0][0],  9.999257249850045e-1, 1.0e-12, "sla::prec", "00", status);
    vvd(mat[0][1], -1.117719859160180e-2, 1.0e-12, "sla::prec", "01", status);
    vvd(mat[0][2], -4.859500474027002e-3, 1.0e-12, "sla::prec", "02", status);
    vvd(mat[1][0],  1.117719858025860e-2, 1.0e-12, "sla::prec", "10", status);
    vvd(mat[1][1],  9.999375327960091e-1, 1.0e-12, "sla::prec", "11", status);
    vvd(mat[1][2], -2.716114374174549e-5, 1.0e-12, "sla::prec", "12", status);
    vvd(mat[2][0],  4.859500500117173e-3, 1.0e-12, "sla::prec", "20", status);
    vvd(mat[2][1], -2.715647545167383e-5, 1.0e-12, "sla::prec", "21", status);
    vvd(mat[2][2],  9.999881921889954e-1, 1.0e-12, "sla::prec", "22", status);

    precl(1925.0, 1975.0, mat);
    vvd(mat[0][0],  9.999257331781050e-1, 1.0e-12, "sla::precl", "00", status);
    vvd(mat[0][1], -1.117658038434041e-2, 1.0e-12, "sla::precl", "01", status);
    vvd(mat[0][2], -4.859236477249598e-3, 1.0e-12, "sla::precl", "02", status);
    vvd(mat[1][0],  1.117658037299592e-2, 1.0e-12, "sla::precl", "10", status);
    vvd(mat[1][1],  9.999375397061558e-1, 1.0e-12, "sla::precl", "11", status);
    vvd(mat[1][2], -2.715816653174189e-5, 1.0e-12, "sla::precl", "12", status);
    vvd(mat[2][0],  4.859236503342703e-3, 1.0e-12, "sla::precl", "20", status);
    vvd(mat[2][1], -2.715349745834860e-5, 1.0e-12, "sla::precl", "21", status);
    vvd(mat[2][2],  9.999881934719490e-1, 1.0e-12, "sla::precl", "22", status);
}

// tests sla::prenut() function
static void t_prenut(bool& status) {
    Matrix<double> mat;
    prenut(1985.0, 50123.4567, mat);
    vvd(mat[0][0],  9.999962358680738e-1, 1.0e-12, "sla::prenut", "00", status);
    vvd(mat[0][1], -2.516417057665452e-3, 1.0e-12, "sla::prenut", "01", status);
    vvd(mat[0][2], -1.093569785342370e-3, 1.0e-12, "sla::prenut", "02", status);
    vvd(mat[1][0],  2.516462370370876e-3, 1.0e-12, "sla::prenut", "10", status);
    vvd(mat[1][1],  9.999968329010883e-1, 1.0e-12, "sla::prenut", "11", status);
    vvd(mat[1][2],  4.006159587358310e-5, 1.0e-12, "sla::prenut", "12", status);
    vvd(mat[2][0],  1.093465510215479e-3, 1.0e-12, "sla::prenut", "20", status);
    vvd(mat[2][1], -4.281337229063151e-5, 1.0e-12, "sla::prenut", "21", status);
    vvd(mat[2][2],  9.999994012499173e-1, 1.0e-12, "sla::prenut", "22", status);
}

// tests sla::dsep(), sla::dsepv(), sla::sep(), and sla::sepv() functions
static void t_sep(bool& status) {
    const Vector<float> vf1 = {1.0f, 0.1f, 0.2f};
    const Vector<float> vf2 = {-3.0f, 1.0e-3f, 0.2f};
    const Vector<double> vd1 = {1.0, 0.1, 0.2};
    const Vector<double> vd2 = {-3.0, 1.0e-3, 0.2};

    Spherical<double> sd1, sd2;
    dcc2s (vd1, sd1);
    dcc2s (vd2, sd2);
    const Spherical<float> sf1 = {(float) sd1.get_ra(), (float) sd1.get_dec()};
    const Spherical<float> sf2 = {(float) sd2.get_ra(), (float) sd2.get_dec()};

    vvd(dsep(sd1, sd2), 2.8603919190246608, 1.0e-7, "sla::dsep", "", status);
    vvd(sep(sf1, sf2), 2.8603919190246608, 1.0e-4, "sla::sep", "", status);
    vvd(dsepv(vd1, vd2), 2.8603919190246608, 1.0e-7, "sla::dsepv", "", status);
    vvd(sepv(vf1, vf2), 2.8603919190246608, 1.0e-4, "sla::sepv", "", status);
}

// tests sla::pa() function
static void t_pa(bool& status) {
    vvd(pa({-1.567, 1.5123}, 0.987), -1.486288540423851, 1.0e-12, "sla::pa", "", status);
    vvd(pa({0.0, 0.789}, 0.789), 0.0, 0.0, "sla::pa", "zenith", status);
}

// tests sla::rcc() function
static void t_rcc(bool& status) {
    vvd(rcc(48939.123, 0.76543, 5.0123, 5525.242, 3190.0),
        -1.280131613589158e-3, 1.0e-15, "sla::rcc", "", status);
}

// tests sla::gmst() and sla::gmsta functions
static void t_gmst(bool& status) {
    vvd(gmst(43999.999), 3.9074971356487318, 1.0e-9, "sla::gmst", "", status);
    vvd(gmsta(43999.0, 0.999), 3.9074971356487318, 1.0e-12, "sla::gmsta", "", status);
}

// tests sla::rcc() function
static void t_prebn(bool& status) {
    Matrix<double> mat;
    prebn(1925.0, 1975.0, mat);
    vvd(mat[0][0],  9.999257613786738e-1, 1.0e-12, "sla::prebn", "00", status);
    vvd(mat[0][1], -1.117444640880939e-2, 1.0e-12, "sla::prebn", "01", status);
    vvd(mat[0][2], -4.858341150654265e-3, 1.0e-12, "sla::prebn", "02", status);
    vvd(mat[1][0],  1.117444639746558e-2, 1.0e-12, "sla::prebn", "10", status);
    vvd(mat[1][1],  9.999375635561940e-1, 1.0e-12, "sla::prebn", "11", status);
    vvd(mat[1][2], -2.714797892626396e-5, 1.0e-12, "sla::prebn", "12", status);
    vvd(mat[2][0],  4.858341176745641e-3, 1.0e-12, "sla::prebn", "20", status);
    vvd(mat[2][1], -2.714330927085065e-5, 1.0e-12, "sla::prebn", "21", status);
    vvd(mat[2][2],  9.999881978224798e-1, 1.0e-12, "sla::prebn", "22", status);
}

// tests sla::rcc() function
static void t_preces(bool& status) {
    Spherical<double> pos = {6.28, -1.123};
    preces(CAT_FK4, 1925.0, 1950.0, pos);
    vvd(pos.get_ra(),  0.002403604864728447, 1.0e-12, "sla::preces", "ra", status);
    vvd(pos.get_dec(), -1.120570643322045, 1.0e-12, "sla::preces", "dec", status);

    pos.set_ra(0.0123);
    pos.set_dec(1.0987);
    preces(CAT_FK5, 2050.0, 1990.0, pos);
    vvd(pos.get_ra(), 6.282003602708382, 1.0e-12, "sla::preces", "ra", status);
    vvd(pos.get_dec(), 1.092870326188383, 1.0e-12, "sla::preces", "dec", status);
}

// tests sla::rcc() function
static void t_supgal(bool& status) {
    const Spherical<double> sgal = {6.1, -1.4};
    Spherical<double> gal;

    supgal(sgal, gal);
    vvd(gal.get_longitude(), 3.798775860769474, 1.0e-12, "sla::supgal", "longitude", status);
    vvd(gal.get_latitude(), -0.1397070490669407, 1.0e-12, "sla::supgal", "latitude", status);
}

// tests sla::rverot(), sla::rvgalc(), sla::rvlg(), sla::rvlsrd(), and sla::rvlsrk() functions
static void t_rv(bool& status) {
    vvd(rverot(-0.777, {5.67, -0.3}, 3.19), -0.1948098355075913, 1.0e-6, "sla::rverot", "", status);
    vvd(rvgalc({1.11, -0.99}), 158.9630759840254, 1.0e-3, "sla::rvgalc", "", status);
    vvd(rvlg({3.97, 1.09}), -197.818762175363, 1.0e-3, "sla::rvlg", "", status);
    vvd(rvlsrd({6.01, 0.1}), -4.082811335150567, 1.0e-4, "sla::rvlsrd", "", status);
    vvd(rvlsrk({6.01, 0.1}), -5.925180579830265, 1.0e-4, "sla::rvlsrk", "", status);
}

// tests sla::cc62s() and sla::dc62s() functions
static void t_cc62s(bool& status) {
    const VectorPV<float> v({100.0f, -50.0f, 25.0f}, {-0.1f, 0.2f, 0.7f});
    SphericalPV<float> s;
    cc62s(v, s);
    vvd(s.get_longitude(), -0.4636476090008061, 1.0e-6, "sla::cc62s", "longitude", status);
    vvd(s.get_latitude(), 0.2199879773954594, 1.0e-6, "sla::cc62s", "latitude", status);
    vvd(s.get_dist(), 114.564392373896, 1.0e-3, "sla::cc62s", "dist", status);
    vvd(s.get_dlong(), 0.001200000000000000, 1.0e-9, "sla::cc62s", "dlong", status);
    vvd(s.get_dlat(), 0.006303582107999407, 1.0e-8, "sla::cc62s", "dlat", status);
    vvd(s.get_ddist(), -0.02182178902359925, 1.0e-7, "sla::cc62s", "ddist", status);

    const VectorPV<double> dv({100.0, -50.0, 25.0}, {-0.1, 0.2, 0.7});
    SphericalPV<double> ds;
    dc62s(dv, ds);
    vvd(ds.get_longitude(), -0.4636476090008061, 1.0e-6, "sla::dc62s", "longitude", status);
    vvd(ds.get_latitude(), 0.2199879773954594, 1.0e-6, "sla::dc62s", "latitude", status);
    vvd(ds.get_dist(), 114.564392373896, 1.0e-9, "sla::dc62s", "dist", status);
    vvd(ds.get_dlong(), 0.001200000000000000, 1.0e-15, "sla::dc62s", "dlong", status);
    vvd(ds.get_dlat(), 0.006303582107999407, 1.0e-14, "sla::dc62s", "dlat", status);
    vvd(ds.get_ddist(), -0.02182178902359925, 1.0e-13, "sla::dc62s", "ddist", status);
}

// tests sla::cs2c6() and sla::ds2c6() functions
static void t_cs2c6(bool& status) {
    VectorPV<float> spv;
    cs2c6({{{-3.21f, 0.123f}, 0.456f}, {{-7.8e-6f, 9.01e-6f}, -1.23e-5f}}, spv);
    vvd(spv.get_x(), -0.4514964673880165, 1.0e-6, "sla::cs2c6", "x", status);
    vvd(spv.get_y(), 0.03093394277342585, 1.0e-6, "sla::cs2c6", "y", status);
    vvd(spv.get_z(), 0.05594668105108779, 1.0e-6, "sla::cs2c6", "z", status);
    vvd(spv.get_dx(), 1.292270850663260e-5, 1.0e-6, "sla::cs2c6", "xd", status);
    vvd(spv.get_dy(), 2.652814182060692e-6, 1.0e-6, "sla::cs2c6", "yd", status);
    vvd(spv.get_dz(), 2.568431853930293e-6, 1.0e-6, "sla::cs2c6", "zd", status);

    VectorPV<double> dpv;
    ds2c6({{{-3.21, 0.123}, 0.456}, {{-7.8e-6, 9.01e-6}, -1.23e-5}}, dpv);
    vvd(dpv.get_x(), -0.4514964673880165, 1.0e-12, "sla::ds2c6", "x", status);
    vvd(dpv.get_y(), 0.03093394277342585, 1.0e-12, "sla::ds2c6", "y", status);
    vvd(dpv.get_z(), 0.05594668105108779, 1.0e-12, "sla::ds2c6", "z", status);
    vvd(dpv.get_dx(), 1.292270850663260e-5, 1.0e-12, "sla::ds2c6", "xd", status);
    vvd(dpv.get_dy(), 2.652814182060692e-6, 1.0e-12, "sla::ds2c6", "yd", status);
    vvd(dpv.get_dz(), 2.568431853930293e-6, 1.0e-12, "sla::ds2c6", "zd", status);
}

// tests sla::etrms() function
static void t_etrms(bool& status) {
    Vector<double> et;
    etrms (1976.9, et);
    vvd(et[0], -1.621617102537041e-6, 1.0e-18, "sla::etrms", "x", status);
    vvd(et[1], -3.310070088507914e-7, 1.0e-18, "sla::etrms", "y", status);
    vvd(et[2], -1.435296627515719e-7, 1.0e-18, "sla::etrms", "z", status);
}

// tests sla::addet() and sla::subet() functions
static void t_addet(bool& status) {
    const Spherical<double> dir {2.0, -1.0};
    Spherical<double> dir1, dir2;
    constexpr double be = 1975.0;

    addet(dir, be, dir1);
    vvd(dir1.get_ra() - dir.get_ra(), 2.983864874295250e-6, 1.0e-12, "sla::addet", "ra", status);
    vvd(dir1.get_dec() - dir.get_dec(), 2.379650804185118e-7, 1.0e-12, "sla::addet", "dec", status);

    subet(dir1, be, dir2);
    vvd(dir2.get_ra() - dir.get_ra(), 0.0, 1.0e-12, "sla::subet", "ra", status);
    vvd(dir2.get_dec() - dir.get_dec(), 0.0, 1.0e-12, "sla::subet", "dec", status);
}

// tests sla::pvobs() and (indirectly) sla::geoc() functions
static void t_pvobs(bool& status) {
    VectorPV<double> pv;

    pvobs(0.5123, 3001.0, -0.567, pv);
    vvd(pv.get_x(), 0.3138647803054939e-4, 1.0e-16,  "sla::pvobs", "x", status);
    vvd(pv.get_y(),-0.1998515596527082e-4, 1.0e-16,  "sla::pvobs", "y", status);
    vvd(pv.get_z(), 0.2078572043443275e-4, 1.0e-16,  "sla::pvobs", "z", status);
    vvd(pv.get_dx(), 0.1457340726851264e-8, 1.0e-20, "sla::pvobs", "dx", status);
    vvd(pv.get_dy(), 0.2288738340888011e-8, 1.0e-20, "sla::pvobs", "dy", status);
    vvd(pv.get_dz(), 0.0, 0.0, "sla::pvobs", "dz", status);
}

// tests sla::pcd() and sla::unpcd() functions
static void t_pcd(bool& status) {
    constexpr double disco = 178.585;
    double x = 0.0123;
    double y = -0.00987;

    pcd(disco, x, y);
    vvd(x, 0.01284630845735895, 1.0e-14, "sla::pcd", "x", status);
    vvd(y, -0.01030837922553926, 1.0e-14, "sla::pcd", "y", status);

    unpcd(disco, x, y);
    vvd(x, 0.0123, 1.0e-14, "sla::unpcd", "x", status);
    vvd(y, -0.00987, 1.0e-14, "sla::unpcd", "y", status);
}

// tests sla::eqeqx() function
static void t_eqeqx(bool& status) {
    vvd(eqeqx(41234.0), 5.376047445838358596e-5, 1.0e-17, "sla::eqeqx", "", status);
}

// tests sla::eqecl() function
static void t_eqecl(bool& status) {
    Spherical<double> dir;
    eqecl({0.789, -0.123}, 46555.0, dir);

    vvd(dir.get_longitude(), 0.7036566430349022, 1.0e-12, "sla::eqecl", "longitude", status);
    vvd(dir.get_latitude(), -0.4036047164116848, 1.0e-12, "sla::eqecl", "latitude", status);
}

// tests sla::eqgal() function
static void t_eqgal(bool& status) {
    Spherical<double> gal;
    eqgal({5.67, -1.23}, gal);

    vvd(gal.get_longitude(), 5.612270780904526, 1.0e-12, "sla::eqgal", "longitude", status);
    vvd(gal.get_latitude(), -0.6800521449061520, 1.0e-12, "sla::eqgal", "latitude", status);
}

// tests sla::galeq() function
static void t_galeq(bool& status) {
    Spherical<double> dir;
    galeq({5.67, -1.23}, dir);

    vvd(dir.get_longitude(), 0.04729270418071426, 1.0e-12, "sla::galeq", "ra", status);
    vvd(dir.get_latitude(), -0.7834003666745548, 1.0e-12, "sla::galeq", "dec", status);
}

// tests sla::fitxy(), sla::pxy(), sla::invf(), sla::xy2xy(), and sla::dcmpf() functions
static void t_fitxy(bool& status) {
    constexpr int NUM_POINTS = 8;
    static const XYSamples expected = {
        {-23.4, -12.1}, {32.0, -15.3}, {10.9, 23.7}, {-3.0, 16.1},
        {45.0, 32.5}, {8.6, -17.0}, {15.3, 10.0}, {121.7, -3.8}
    };
    static const XYSamples measured = {
        {-23.41, 12.12}, {32.03, 15.34}, {10.93, -23.72}, {-3.01, -16.10},
        {44.90, -32.46}, {8.55, 17.02}, {15.31, -10.07}, {120.92, 3.81}
    };
    double predicted[NUM_POINTS][2], x_rms, y_rms, rms, x2, y2, xz, yz, xs, ys, perp, orient;
    FitCoeffs model, inverse;

    // fit a 4-coeff linear model to relate two sets of (x,y) coordinates
    FITStatus result = fitxy(true, NUM_POINTS, expected, measured, model);
    vvd(model[0], -7.938263381515947e-3, 1.0e-12, "sla::fitxy", "4/0", status);
    vvd(model[1], 1.004640925187200, 1.0e-12, "sla::fitxy", "4/1", status);
    vvd(model[2], 3.976948048238268e-4, 1.0e-12, "sla::fitxy", "4/2", status);
    vvd(model[3], -2.501031681585021e-2, 1.0e-12, "sla::fitxy", "4/3", status);
    vvd(model[4], 3.976948048238268e-4, 1.0e-12, "sla::fitxy", "4/4", status);
    vvd(model[5], -1.004640925187200, 1.0e-12, "sla::fitxy", "4/5", status);
    viv(result, FIT_OK, "sla::fitxy", "4/result", status);

    // same but 6-coeff
    result = fitxy(false, NUM_POINTS, expected, measured, model);
    vvd(model[0], -2.617232551841476e-2, 1.0e-12, "sla::fitxy", "6/0", status);
    vvd(model[1], 1.005634905041421, 1.0e-12, "sla::fitxy", "6/1", status);
    vvd(model[2], 2.133045023329208e-3, 1.0e-12, "sla::fitxy", "6/2", status);
    vvd(model[3], 3.846993364417779909e-3, 1.0e-12, "sla::fitxy", "6/3", status);
    vvd(model[4], 1.301671386431460e-4, 1.0e-12, "sla::fitxy", "6/4", status);
    vvd(model[5], -0.9994827065693964, 1.0e-12, "sla::fitxy", "6/5", status);
    viv(result, FIT_OK, "sla::fitxy", "6/result", status);

    // compute predicted coordinates and residuals
    pxy(NUM_POINTS, expected, measured, model, predicted, x_rms, y_rms, rms);
    vvd(predicted[0][0], -23.542232946855340, 1.0e-12, "sla::pxy", "x0", status);
    vvd(predicted[0][1], -12.11293062297230597, 1.0e-12, "sla::pxy", "y0", status);
    vvd(predicted[1][0], 32.217034593616180, 1.0e-12, "sla::pxy", "x1", status);
    vvd(predicted[1][1], -15.324048471959370, 1.0e-12, "sla::pxy", "y1", status);
    vvd(predicted[2][0], 10.914821358630950, 1.0e-12, "sla::pxy", "x2", status);
    vvd(predicted[2][1], 23.712999520015880, 1.0e-12, "sla::pxy", "y2", status);
    vvd(predicted[3][0], -3.087475414568693, 1.0e-12, "sla::pxy", "x3", status);
    vvd(predicted[3][1], 16.09512676604438414, 1.0e-12, "sla::pxy", "y3", status);
    vvd(predicted[4][0], 45.05759626938414666, 1.0e-12, "sla::pxy", "x4", status);
    vvd(predicted[4][1], 32.45290015313210889, 1.0e-12, "sla::pxy", "y4", status);
    vvd(predicted[5][0], 8.608310538882801, 1.0e-12, "sla::pxy", "x5", status);
    vvd(predicted[5][1], -17.006235743411300, 1.0e-12, "sla::pxy", "y5", status);
    vvd(predicted[6][0], 15.348618307280820, 1.0e-12, "sla::pxy", "x6", status);
    vvd(predicted[6][1], 10.07063070741086835, 1.0e-12, "sla::pxy", "y6", status);
    vvd(predicted[7][0], 121.5833272936291482, 1.0e-12, "sla::pxy", "x7", status);
    vvd(predicted[7][1], -3.788442308260240, 1.0e-12, "sla::pxy", "y7", status);
    vvd(x_rms ,0.1087247110488075, 1.0e-13, "sla::pxy", "x_rms", status);
    vvd(y_rms, 0.03224481175794666, 1.0e-13, "sla::pxy", "y_rms", status);
    vvd(rms, 0.1134054261398109, 1.0e-13, "sla::pxy", "rms", status);

    // invert the model
    bool success = invf(model, inverse);
    vvd(inverse[0], 0.02601750208015891, 1.0e-12, "sla::invf", "0", status);
    vvd(inverse[1], 0.9943963945040283, 1.0e-12, "sla::invf", "1", status);
    vvd(inverse[2], 0.002122190075497872, 1.0e-12, "sla::invf", "2", status);
    vvd(inverse[3], 0.003852372795357474353, 1.0e-12, "sla::invf", "3", status);
    vvd(inverse[4], 0.0001295047252932767, 1.0e-12, "sla::invf", "4", status);
    vvd(inverse[5], -1.000517284779212, 1.0e-12, "sla::invf", "5", status);
    viv(success, (int) true, "sla::invf", "success", status);

    // transform one [x,y]
    xy2xy(44.5, 32.5, model, x2, y2);
    vvd(x2, 44.793904912083030, 1.0e-11, "sla::xy2xy", "x", status);
    vvd(y2, -32.473548532471330, 1.0e-11, "sla::xy2xy", "y", status);

    // decompose the fit into scales etc.
    dcmpf(model, xz, yz, xs, ys, perp, orient);
    vvd(xz, -0.0260175020801628646, 1.0e-12, "sla::dcmpf", "xz", status);
    vvd(yz, -0.003852372795357474353, 1.0e-12, "sla::dcmpf", "yz", status);
    vvd(xs, -1.00563491346569, 1.0e-12, "sla::dcmpf", "xs", status);
    vvd(ys, 0.999484982684761, 1.0e-12, "sla::dcmpf", "ys", status);
    vvd(perp,-0.002004707996156263, 1.0e-12, "sla::dcmpf", "perp", status);
    vvd(orient, 3.14046086182333, 1.0e-12, "sla::dcmpf", "orient", status);
}

// tests sla::pm() function
static void t_pm(bool& status) {
    Spherical<double> dir;
    pm({5.43, -0.87}, {-0.33e-5, 0.77e-5},
        0.7, 50.3*365.2422/365.25, 1899.0, 1943.0, dir);
    vvd(dir.get_ra(), 5.429855087793875, 1.0e-12, "sla::pm", "ra", status);
    vvd(dir.get_dec(), -0.8696617307805072, 1.0e-12, "sla::pm", "dec", status);
}

// tests sla::earth() function
static void t_earth(bool& status) {
    VectorPV<float> pv;
    earth(1978, 174, 0.87f, pv);

    vvd(pv.get_x(), 3.590867086e-2, 1.0e-6, "sla::earth", "x", status);
    vvd(pv.get_y(), -9.319285116e-1, 1.0e-6, "sla::earth", "y", status);
    vvd(pv.get_z(), -4.041039435e-1, 1.0e-6, "sla::earth", "z", status);
    vvd(pv.get_dx(), 1.956930055e-7, 1.0e-13, "sla::earth", "dx", status);
    vvd(pv.get_dy(), 5.743797400e-9, 1.0e-13, "sla::earth", "dy", status);
    vvd(pv.get_dz(), 2.512001677e-9, 1.0e-13, "sla::earth", "dz", status);
}

// tests sla::ecor() function
static void t_ecor(bool& status) {
    float velocity, lt;
    ecor({2.345f, -0.567f}, 1995, 306, 0.037f, velocity, lt);

    vvd(velocity, -19.182460, 1.0e-3, "sla::ecor", "velocity", status);
    vvd(lt, -120.36632, 1.0e-2, "sla::ecor", "lt", status);
}

// tests sla::ecleq() function
static void t_ecleq(bool& status) {
    Spherical<double> dir;
    ecleq({1.234, -0.123}, 43210.0, dir);

    vvd(dir.get_ra(), 1.229910118208851, 1.0e-12, "sla::ecleq", "ra", status);
    vvd(dir.get_dec(), 0.2638461400411088, 1.0e-12, "sla::ecleq", "dec", status);
}

// tests sla::polmo() function
static void t_polmo(bool& status) {
    double t_long, t_phi, d_az;
    polmo(0.7, -0.5, 1.0e-6, -2.0e-6, t_long, t_phi, d_az);

    vvd(t_long,  0.7000004837322044, 1.0e-12, "sla::polmo", "long", status);
    vvd(t_phi, -0.4999979467222241, 1.0e-12, "sla::polmo", "phi", status);
    vvd(d_az,  1.008982781275728e-6, 1.0e-12, "sla::polmo", "az", status);
}

// tests sla::galsup() function
static void t_galsup(bool& status) {
    Spherical<double> supergalactic;
    galsup({6.1, -1.4}, supergalactic);

    vvd(supergalactic.get_longitude(), 4.567933268859171, 1.0e-12, "sla::galsup", "long", status);
    vvd(supergalactic.get_latitude(), -0.01862369899731829, 1.0e-12, "sla::galsup", "lat", status);
}

// tests sla::s2tp(), sla::ds2tp(), sla::dtps2c(), sla::tp2s(), sla::dtp2s(), and sla::tps2c() functions
static void t_tp(bool& status) {
    const float r0 = 3.1f;
    const float d0 = -0.9f;
    const float r1 = r0 + 0.2f;
    const float d1 = d0 - 0.1f;
    float x, y;
    TPPStatus result = s2tp({r1, d1}, {r0, d0}, x, y);
    vvd(x, 0.1086112301590404, 1.0e-6, "sla::s2tp", "x", status);
    vvd(y, -0.1095506200711452, 1.0e-6, "sla::s2tp", "y", status);
    viv(result, TPP_OK, "sla::s2tp", "result", status);
    Spherical<float> point;
    tp2s(x, y, {r0, d0}, point);
    vvd(point.get_ra() - r1, 0.0, 1.0e-6, "sla::tp2s", "ra", status);
    vvd(point.get_dec() - d1, 0.0, 1.0e-6, "sla::tp2s", "dec", status);
    Spherical<float> s1, s2;
    int n = tps2c(x, y, point, s1, s2);
    vvd(s1.get_ra(),  3.1, 1.0e-6, "sla::tps2c", "ra1", status);
    vvd(s1.get_dec(), -0.9, 1.0e-6, "sla::tps2c", "dec1", status);
    vvd(s2.get_ra(), 0.3584073464102072, 1.0e-6, "sla::tps2c", "ra2", status);
    vvd(s2.get_dec(), -2.023361658234722, 1.0e-6, "sla::tps2c", "dec2", status);
    viv(n, 1, "sla::tps2c", "n", status);

    const double dr0 = 3.1;
    const double dd0 = -0.9;
    const double dr1 = dr0 + 0.2;
    const double dd1 = dd0 - 0.1;
    double dx, dy;
    result = ds2tp({dr1, dd1}, {dr0, dd0}, dx, dy);
    vvd(dx, 0.1086112301590404, 1.0e-12, "sla::ds2tp", "x", status);
    vvd(dy, -0.1095506200711452, 1.0e-12, "sla::ds2tp", "y", status);
    viv(result, TPP_OK, "sla::ds2tp", "result", status);
    Spherical<double> dpoint;
    dtp2s(dx, dy, {dr0, dd0}, dpoint);
    vvd(dpoint.get_ra() - dr1, 0.0, 1.0e-12, "sla::dtp2s", "ra", status);
    vvd(dpoint.get_dec() - dd1, 0.0, 1.0e-12, "sla::dtp2s", "dec", status);
    Spherical<double> ds1, ds2;
    n = dtps2c(dx, dy, dpoint, ds1, ds2);
    vvd(ds1.get_ra(),  3.1, 1.0e-12, "sla::dtps2c", "ra1", status);
    vvd(ds1.get_dec(), -0.9, 1.0e-12, "sla::dtps2c", "dec1", status);
    vvd(ds2.get_ra(), 0.3584073464102072, 1.0e-12, "sla::dtps2c", "ra2", status);
    vvd(ds2.get_dec(), -2.023361658234722, 1.0e-12, "sla::dtps2c", "dec2", status);
    viv(n, 1, "sla::dtps2c", "n", status);
}

// tests sla::tp2v(), sla::v2tp(), sla::tpv2c(), sla::dtp2v(), sla::dv2tp(), and sla::dtpv2c() functions
static void t_tpv(bool& status) {
    float fr_xi, fr_eta;
    Vector<float> fr_v, fs_v1, fs_v2;
    double dr_xi, dr_eta;
    Vector<double> dr_v, ds_v1, ds_v2;

    const double d_xi = -0.1;
    const double d_eta = 0.055;
    const float f_xi = (float) d_xi;
    const float f_eta = (float) d_eta;

    double x = -0.7;
    double y = -0.13;
    double z = std::sqrt(1.0 - x * x - y * y);
    const Vector<float> f_v = {(float) x, (float) y, (float) z};
    const Vector<double> d_v = {x, y, z};

    x = -0.72;
    y = -0.16;
    z = std::sqrt(1.0 - x * x - y * y);
    const Vector<float> f_v0 = {(float) x, (float) y, (float) z};
    const Vector<double> d_v0 = {x, y, z};

    tp2v(f_xi, f_eta, f_v0, fr_v);
    vvd(fr_v[0], -0.700887428128, 1.0e-6, "sla::tp2v", "v0", status);
    vvd(fr_v[1], -0.05397407, 1.0e-6, "sla::tp2v", "v1", status);
    vvd(fr_v[2], 0.711226836562, 1.0e-6, "sla::tp2v", "v2", status);

    dtp2v(d_xi, d_eta, d_v0, dr_v);
    vvd(dr_v[0], -0.7008874281280771, 1.0e-13, "sla::dtp2v", "v0", status);
    vvd(dr_v[1], -0.05397406827952735, 1.0e-13, "sla::dtp2v", "v1", status);
    vvd(dr_v[2], 0.7112268365615617, 1.0e-13, "sla::dtp2v", "v2", status);

    TPPStatus result = v2tp(f_v, f_v0, fr_xi, fr_eta);
    vvd(fr_xi, -0.02497229197, 1.0e-6, "sla::v2tp", "d_xi", status);
    vvd(fr_eta, 0.03748140764, 1.0e-6, "sla::v2tp", "d_eta", status);
    viv(result, TPP_OK, "sla::v2tp", "result", status);

    result = dv2tp(d_v, d_v0, dr_xi, dr_eta);
    vvd(dr_xi, -0.02497229197023852, 1.0e-13, "sla::dv2tp", "d_xi", status);
    vvd(dr_eta, 0.03748140764224765, 1.0e-13, "sla::dv2tp", "d_eta", status);
    viv(result, TPP_OK, "sla::dv2tp", "result", status);

    int n = tpv2c(f_xi, f_eta, f_v, fs_v1, fs_v2);
    vvd(fs_v1[0], -0.7074573732537283, 1.0e-6, "sla::tpv2c", "v1:0", status);
    vvd(fs_v1[1], -0.2372965765309941, 1.0e-6, "sla::tpv2c", "v1:1", status);
    vvd(fs_v1[2], 0.6657284730245545, 1.0e-6, "sla::tpv2c", "v1:2", status);
    vvd(fs_v2[0], -0.6680480104758149, 1.0e-6, "sla::tpv2c", "v2:0", status);
    vvd(fs_v2[1], -0.02915588494045333, 1.0e-6, "sla::tpv2c", "v2:1", status);
    vvd(fs_v2[2], 0.7435467638774610, 1.0e-6, "sla::tpv2c", "v2:2", status);
    viv(n, 1, "sla::tpv2c", "n", status);

    n = dtpv2c(d_xi, d_eta, d_v, ds_v1, ds_v2);
    vvd(ds_v1[0], -0.7074573732537283, 1.0e-13, "sla::dtpv2c", "v1:0", status);
    vvd(ds_v1[1], -0.2372965765309941, 1.0e-13, "sla::dtpv2c", "v1:1", status);
    vvd(ds_v1[2], 0.6657284730245545, 1.0e-13, "sla::dtpv2c", "v1:2", status);
    vvd(ds_v2[0], -0.6680480104758149, 1.0e-13, "sla::dtpv2c", "v2:0", status);
    vvd(ds_v2[1], -0.02915588494045333, 1.0e-13, "sla::dtpv2c", "v2:1", status);
    vvd(ds_v2[2], 0.7435467638774610, 1.0e-13, "sla::dtpv2c", "v2:2", status);
    viv(n, 1, "sla::dtpv2c", "n", status);
}

// tests sla::combn() and sla::permut() functions
static void t_percom(bool& status) {
    CPStatus result;
    int list[3] = {0};
    for (int i = 0; i < 11; i++) {
        result = combn(3, 5, list);
    }
    viv(result, CPS_NO_MORE, "sla::combn", "result", status);
    viv(list[0], 1, "sla::combn", "list:0", status);
    viv(list[1], 2, "sla::combn", "list:1", status);
    viv(list[2], 3, "sla::combn", "list:2", status);

    int state[4] = {-1}, order[4];
    for (int j = 0; j < 25; j++) {
        result = permut(4, state, order);
    }
    viv(result, CPS_NO_MORE, "sla::permut", "result", status);
    viv(order[0], 4, "sla::permut", "order:0", status);
    viv(order[1], 3, "sla::permut", "order:1", status);
    viv(order[2], 2, "sla::permut", "order:2", status);
    viv(order[3], 1, "sla::permut", "order:3", status);
}

// tests sla::evp() and sla::epv() functions
static void t_evp(bool& status) {
    Vector<double> bvelo, bpos, hvelo, hpos;

    evp(50100.0, 1990.0, bvelo, bpos, hvelo, hpos);
    vvd(bvelo[0], -1.807210068604058436e-7, 1e-14, "sla::evp", "bvelo:x", status);
    vvd(bvelo[1], -8.385891022440320e-8, 1e-14, "sla::evp", "bvelo:y", status);
    vvd(bvelo[2], -3.635846882638055e-8, 1e-14, "sla::evp", "bvelo:z", status);
    vvd(bpos[0], -0.4515615297360333, 1e-7, "sla::evp", "bpos:x", status);
    vvd(bpos[1],  0.8103788166239596, 1e-7, "sla::evp", "bpos:y", status);
    vvd(bpos[2],  0.3514505204144827, 1e-7, "sla::evp", "bpos:z", status);
    vvd(hvelo[0], -1.806354061156890855e-7, 1e-14, "sla::evp", "hvelo:x", status);
    vvd(hvelo[1], -8.383798678086174e-8, 1e-14, "sla::evp", "hvelo:y", status);
    vvd(hvelo[2], -3.635185843644782e-8, 1e-14, "sla::evp", "hvelo:z", status);
    vvd(hpos[0], -0.4478571659918565, 1e-7, "sla::evp", "hpos:x", status);
    vvd(hpos[1],  0.8036439916076232, 1e-7, "sla::evp", "hpos:y", status);
    vvd(hpos[2],  0.3484298459102053, 1e-7, "sla::evp", "hpos:z", status);

    epv(53411.52501161, hpos, hvelo, bpos, bvelo);
    vvd(hpos[0], -0.7757238809297653, 1.0e-12, "sla::epv", "hpos:x", status);
    vvd(hpos[1], +0.5598052241363390, 1.0e-12, "sla::epv", "hpos:y", status);
    vvd(hpos[2], +0.2426998466481708, 1.0e-12, "sla::epv", "hpos:z", status);
    vvd(hvelo[0], -0.0109189182414732, 1.0e-12, "sla::epv", "hvelo:x", status);
    vvd(hvelo[1], -0.0124718726844084, 1.0e-12, "sla::epv", "hvelo:y", status);
    vvd(hvelo[2], -0.0054075694180650, 1.0e-12, "sla::epv", "hvelo:z", status);
    vvd(bpos[0], -0.7714104440491060, 1.0e-12, "sla::epv", "bpos:x", status);
    vvd(bpos[1], +0.5598412061824225, 1.0e-12, "sla::epv", "bpos:y", status);
    vvd(bpos[2], +0.2425996277722475, 1.0e-12, "sla::epv", "bpos:z", status);
    vvd(bvelo[0], -0.0109187426811683, 1.0e-12, "sla::epv", "bvelo:x", status);
    vvd(bvelo[1], -0.0124652546173285, 1.0e-12, "sla::epv", "bvelo:y", status);
    vvd(bvelo[2], -0.0054047731809662, 1.0e-12, "sla::epv", "bvelo:z", status);
}

// tests sla::eg50() function
static void t_eg50(bool& status) {
    Spherical<double> gal;
    eg50({3.012, 1.234}, gal);
    vvd(gal.get_longitude(), 2.305557953813397, 1.0e-12, "sla::eg50", "l", status);
    vvd(gal.get_latitude(), 0.7903600886585871, 1.0e-12, "sla::eg50", "b", status);
}

// tests sla::ge50() function
static void t_ge50(bool& status) {
    Spherical<double> loc;
    ge50({6.1, -1.55}, loc);
    vvd(loc.get_ra(), 0.1966825219934508, 1.0e-12, "sla::ge50", "ra", status);
    vvd(loc.get_dec(), -0.4924752701678960, 1.0e-12, "sla::ge50", "dec", status);
}

// tests sla::pdq2h() function
static void t_pdq2h(bool& status) {
    bool ha1_valid, ha2_valid;
    double ha1, ha2;

    pdq2h(0.9, 0.2, 0.1, ha1, ha1_valid, ha2, ha2_valid);
    vvd(ha1, 0.1042809894435257, 1.0e-14, "sla::pdq2h", "ha1", status);
    viv(ha1_valid, true, "sla::pdq2h", "ha1_v", status);
    vvd(ha2, 2.997450098818439, 1.0e-13, "sla::pdq2h", "ha2", status);
    viv(ha2_valid, true, "sla::pdq2h", "ha2_v", status);
}

// tests sla::pda2h() function
static void t_pda2h(bool& status) {
    double ha1, ha2;
    bool ha1_valid, ha2_valid;

    pda2h(-0.51, -1.31, 3.1, ha1, ha1_valid, ha2, ha2_valid);
    vvd(ha1, -0.1161784556585304927, 1.0e-14, "sla::pda2h", "ha1", status);
    viv(ha1_valid, true, "sla::pda2h", "ha1_v", status);
    vvd(ha2, -2.984787179226459, 1.0e-13, "sla::pda2h", "ha2", status);
    viv(ha2_valid, true, "sla::pda2h", "ha2_v", status);
}

// tests sla::moon() and sla::dmoon() functions
static void t_moon(bool& status) {
    VectorPV<float> pv;
    moon(1999, 365, 0.9f, pv);
    vvd(pv.get_x(), -2.155729505970773e-3, 1.0e-6, "sla::moon", "x", status);
    vvd(pv.get_y(), -1.538107758633427e-3, 1.0e-6, "sla::moon", "y", status);
    vvd(pv.get_z(), -4.003940552689305e-4, 1.0e-6, "sla::moon", "z", status);
    vvd(pv.get_dx(), 3.629209419071314e-9, 1.0e-12, "sla::moon", "dx", status);
    vvd(pv.get_dy(), -4.989667166259157e-9, 1.0e-12, "sla::moon", "dy", status);
    vvd(pv.get_dz(), -2.160752457288307e-9, 1.0e-12, "sla::moon", "dz", status);

    /*
     * The original FORTRAN implementation of the `T_MOON` subroutine was missing test for the `sla_DMOON`
     * subroutine, even though it was mentioned in its comment ("Test sla_MOON and sla_DMOON routines.").
     *
     * The below test was put together by simply converting the date (365-th day of the year 1999) used in the
     * `sla::moon()` test above to MJD, adding 0.9 `fraction` to it, and feeding to the `sla::dmoon()`. The position
     * part of the result matches that of the `sla::moon()` down to specified accuracy, while tolerance for the
     * velocity part had to be increased from 1e-12 to 1e-11. This does not mean, of course, that the results of
     * `sla::dmoon()` are less precise; quite the opposite is true: below, we're essentially comparing output of the
     * reference FORTRAN implementation of single-precision `sla_MOON` to the actual output of the
     * double-precision `sla::dmoon()` implementing a more complete version of the algorithm.
     */
    VectorPV<double> dpv;
    dmoon(51543.9, dpv);
    vvd(dpv.get_x(), -2.155729505970773e-3, 1.0e-6, "sla::dmoon", "x", status);
    vvd(dpv.get_y(), -1.538107758633427e-3, 1.0e-6, "sla::dmoon", "y", status);
    vvd(dpv.get_z(), -4.003940552689305e-4, 1.0e-6, "sla::dmoon", "z", status);
    vvd(dpv.get_dx(), 3.629209419071314e-9, 1.0e-11 , "sla::dmoon", "dx", status);
    vvd(dpv.get_dy(), -4.989667166259157e-9, 1.0e-11, "sla::dmoon", "dy", status);
    vvd(dpv.get_dz(), -2.160752457288307e-9, 1.0e-11, "sla::dmoon", "dz", status);
}

///////////////////////////////////////////////////////////////////////////////
// MODULE ENTRY POINT
///////////////////////////////////////////////////////////////////////////////

/// Tests all SLALIB functions and procedures.
bool sla_test() {
    bool status = true;
    t_airmas(status);
    t_bear(status);
    t_caf2r(status);
    t_caldj(status);
    t_calyd(status);
    t_djcal(status);
    t_cc2s(status);
    t_cldj(status);
    t_e2h(status);
    t_vecmat(status);
    t_zd(status);
    t_pa(status);
    t_cd2tf(status);
    t_cr2af(status);
    t_cr2tf(status);
    t_ctf2d(status);
    t_ctf2r(status);
    t_dat(status);
    t_range(status);
    t_ranorm(status);
    t_ref(status);
    t_ecmat(status);
    t_dmat(status);
    t_smat(status);
    t_altaz(status);
    t_nut(status);
    t_epj2d(status);
    t_epj(status);
    t_epb2d(status);
    t_epb(status);
    t_epco(status);
    t_prec(status);
    t_prenut(status);
    t_sep(status);
    t_rcc(status);
    t_gmst(status);
    t_prebn(status);
    t_preces(status);
    t_supgal(status);
    t_rv(status);
    t_cc62s(status);
    t_cs2c6(status);
    t_etrms(status);
    t_addet(status);
    t_pvobs(status);
    t_pcd(status);
    t_eqeqx(status);
    t_eqecl(status);
    t_eqgal(status);
    t_galeq(status);
    t_fitxy(status);
    t_pm(status);
    t_earth(status);
    t_ecor(status);
    t_ecleq(status);
    t_polmo(status);
    t_galsup(status);
    t_tp(status);
    t_tpv(status);
    t_percom(status);
    t_evp(status);
    t_eg50(status);
    t_ge50(status);
    t_pdq2h(status);
    t_pda2h(status);
    t_moon(status);
    return status;
}

}
