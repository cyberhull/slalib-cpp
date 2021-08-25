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
    bool singular = dmat(3, mat, vec, det, ws);

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
    double azimuth, az_vel, az_acc, elevation, el_vel, el_acc, p_angle, pa_vel, pa_acc;
    altaz(0.7, -0.7, -0.65,
        azimuth, az_vel, az_acc, elevation, el_vel, el_acc, p_angle, pa_vel, pa_acc);

    vvd(azimuth, 4.400560746660174, 1.0e-12, "sla::altaz", "azimuth", status);
    vvd(az_vel, -0.2015438937145421, 1.0e-13, "sla::altaz", "az_vel", status);
    vvd(az_acc, -0.4381266949668748, 1.0e-13, "sla::altaz", "az_acc", status);
    vvd(elevation, 1.026646506651396, 1.0e-12, "sla::altaz", "elevation", status);
    vvd(el_vel, -0.7576920683826450, 1.0e-13, "sla::altaz", "el_vel", status);
    vvd(el_acc, 0.04922465406857453, 1.0e-14, "sla::altaz", "el_acc", status);
    vvd(p_angle, 1.707639969653937, 1.0e-12, "sla::altaz", "p_angle", status);
    vvd(pa_vel, 0.4717832355365627, 1.0e-13, "sla::altaz", "pa_vel", status);
    vvd(pa_acc, -0.2957914128185515, 1.0e-13, "sla::altaz", "pa_acc", status);
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
    t_prec(status);
    t_prenut(status);
    t_sep(status);
    t_rcc(status);
    t_prebn(status);
    t_preces(status);
    t_supgal(status);
    t_rv(status);
    return status;
}

}
