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
#include <resolv.h>

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

/// Validates a character string result.
static void vcs(const char* str, const char* str_ok, const char* func, const char* test, bool &status) {
    if (std::strcmp(str, str_ok) != 0) {
        err(func, test, status);
        std::printf("\tExpected: '%s'\n", str_ok);
        std::printf("\t  Actual: '%s'\n", str);
    }
}

/// Validates an integer result.
static void viv(const int val, const int val_ok, const char* func, const char* test, bool &status) {
    if (val != val_ok) {
        err(func, test, status);
        std::printf("\tExpected: %d\n", val_ok);
        std::printf("\t  Actual: %d\n", val);
    }
}

/// Validates a long result; in FORTRAN implementation of SLALIB, "long" values are INTEGER*4, so C/C++ `int`s
/// work just fine, there is no need to use `long`s.
static void vlv(const int val, const int val_ok, const char* func, const char* test, bool &status) {
    static_assert(sizeof val >= 4, "`int` values must be at least 32-bit");
    viv(val, val_ok, func, test, status);
}

/// Validates a double precision floating point result.
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

/// Tests sla::airmas() function.
void t_airmas(bool& status) {
    vvd(airmas(1.2354), 3.015698990074724, 1e-12, "sla::airmas", "", status);
}

/**
 * Tests all the 3-component vector and 3x3 matrix procedures:
 *
 *   sla::av2m()   sla::dav2m()
 *   sla::cc2s()   sla::dcc2s()      (NOT tested directly in either FORTRAN or C++ implementations)
 *   sla::cs2c()   sla::dcs2c()      (these two take structures in C++ implementation)
 *   sla::euler()  sla::deuler()
 *   sla::imxv()   sla::dimxv()
 *   sla::m2av()   sla::dm2av()
 *   sla::mxm()    sla::dmxm()
 *   sla::mxv()    sla::dmxv()
 *   sla::vdv()    sla::dvdv()
 *   sla::vn()     sla::dvn()        (these two return values in the C++ implementation)
 *   sla::vxv()    sla::dvxv()
 */
void t_vecmat(bool& status) {
    // tolerances for [most] single and double precision functions' tests, respectively
    constexpr double sp_tolerance = 1.0e-6;
    constexpr double dp_tolerance = 1.0e-12;

    // make a rotation matrix
    vector<float> av;
    matrix<float> rm1;
    av[0] = -0.123f;
    av[1] = 0.0987f;
    av[2] = 0.0654f;
    av2m(av, rm1);
    vvd(double(rm1[0][0]), 0.9930075842721269, sp_tolerance, "sla::av2m", "00", status);
    vvd(double(rm1[0][1]), 0.05902743090199868, sp_tolerance, "sla::av2m", "01", status);
    vvd(double(rm1[0][2]), -0.1022335560329612, sp_tolerance, "sla::av2m", "02", status);
    vvd(double(rm1[1][0]), -0.07113807138648245, sp_tolerance, "sla::av2m", "10", status);
    vvd(double(rm1[1][1]), 0.9903204657727545, sp_tolerance, "sla::av2m", "11", status);
    vvd(double(rm1[1][2]), -0.1191836812279541, sp_tolerance, "sla::av2m", "12", status);
    vvd(double(rm1[2][0]), 0.09420887631983825, sp_tolerance, "sla::av2m", "20", status);
    vvd(double(rm1[2][1]), 0.1256229973879967, sp_tolerance, "sla::av2m", "21", status);
    vvd(double(rm1[2][2]), 0.9875948309655174, sp_tolerance, "sla::av2m", "22", status);

    // make another
    matrix<float> rm2;
    euler("YZY", 2.345E0, -0.333E0, 2.222E0, rm2);
    vvd(double(rm2[0][0]), -0.1681574770810878, sp_tolerance, "sla::euler", "00", status);
    vvd(double(rm2[0][1]), 0.1981362273264315, sp_tolerance, "sla::euler", "01", status);
    vvd(double(rm2[0][2]), 0.9656423242187410, sp_tolerance, "sla::euler", "02", status);
    vvd(double(rm2[1][0]), -0.2285369373983370, sp_tolerance, "sla::euler", "10", status);
    vvd(double(rm2[1][1]), 0.9450659587140423, sp_tolerance, "sla::euler", "11", status);
    vvd(double(rm2[1][2]), -0.2337117924378156, sp_tolerance, "sla::euler", "12", status);
    vvd(double(rm2[2][0]), -0.9589024617479674, sp_tolerance, "sla::euler", "20", status);
    vvd(double(rm2[2][1]), -0.2599853247796050, sp_tolerance, "sla::euler", "21", status);
    vvd(double(rm2[2][2]), -0.1136384607117296, sp_tolerance, "sla::euler", "22", status);

    // combine them
    matrix<float> rm;
    mxm(rm2, rm1, rm);
    vvd(double(rm[0][0]), -0.09010460088585805, sp_tolerance, "sla::mxm", "00", status);
    vvd(double(rm[0][1]), 0.3075993402463796, sp_tolerance, "sla::mxm", "01", status);
    vvd(double(rm[0][2]), 0.9472400998581048, sp_tolerance, "sla::mxm", "02", status);
    vvd(double(rm[1][0]), -0.3161868071070688, sp_tolerance, "sla::mxm", "10", status);
    vvd(double(rm[1][1]), 0.8930686362478707, sp_tolerance, "sla::mxm", "11", status);
    vvd(double(rm[1][2]), -0.3200848543149236, sp_tolerance, "sla::mxm", "12", status);
    vvd(double(rm[2][0]), -0.9444083141897035, sp_tolerance, "sla::mxm", "20", status);
    vvd(double(rm[2][1]), -0.3283459407855694, sp_tolerance, "sla::mxm", "21", status);
    vvd(double(rm[2][2]), 0.01678926022795169, sp_tolerance, "sla::mxm", "22", status);

    // create a vector
    vector<float> v1;
    cs2c({3.0123f, -0.999f}, v1 );
    vvd(double(v1[0]), -0.5366267667260525, sp_tolerance, "sla::cs2c", "X", status);
    vvd(double(v1[1]), 0.06977111097651444, sp_tolerance, "sla::cs2c", "Y", status);
    vvd(double(v1[2]), -0.8409302618566215, sp_tolerance, "sla::cs2c", "Z", status);

    // rotate the vector using the two matrices sequentially
    vector<float> v2, v3;
    mxv(rm1, v1, v2);
    mxv(rm2, v2, v3);
    vvd(double(v3[0]), -0.7267487768696160, sp_tolerance, "sla::mxv", "X", status);
    vvd(double(v3[1]), 0.5011537352639822, sp_tolerance, "sla::mxv", "Y", status);
    vvd(double(v3[2]), 0.4697671220397141, sp_tolerance, "sla::mxv", "Z", status);

    // de-rotate the vector using the combined matrix
    vector<float> v4;
    imxv(rm, v3, v4);
    vvd(double(v4[0]), -0.5366267667260526, sp_tolerance, "sla::imxv", "X", status);
    vvd(double(v4[1]), 0.06977111097651445, sp_tolerance, "sla::imxv", "Y", status);
    vvd(double(v4[2]), -0.8409302618566215, sp_tolerance, "sla::imxv", "Z", status);

    // convert the combined matrix into an axial vector
    vector<float> v5;
    m2av(rm, v5);
    vvd(double(v5[0]), 0.006889040510209034, sp_tolerance, "sla::m2av", "X", status);
    vvd(double(v5[1]), -1.577473205461961, sp_tolerance, "sla::m2av", "Y", status);
    vvd(double(v5[2]), 0.5201843672856759, sp_tolerance, "sla::m2av", "Z", status);

    // multiply the axial vector by a scalar and then normalize
    for (int i = 0; i < 3; i++) {
        v5[i] *= 1000.0f;
    }
    vector<float> v6;
    float vm = vn(v5, v6);
    vvd(double(v6[0]), 0.004147420704640065, sp_tolerance, "sla::vn", "X", status);
    vvd(double(v6[1]), -0.9496888606842218, sp_tolerance, "sla::vn", "Y", status);
    vvd(double(v6[2]), 0.3131674740355448, sp_tolerance, "sla::vn", "Z", status);
    vvd(double(vm), 1661.042127339937, 1.0e-3, "sla::VN", "m", status);

    // calculate dot product with the original vector
    vvd(double(vdv(v6, v1)), -0.3318384698006295, sp_tolerance, "sla::VN", " ", status);

    // calculate cross product with the original vector
    vector<float> v7;
    vxv(v6, v1, v7);
    vvd(double(v7[0]), 0.7767720597123304, sp_tolerance, "sla::vxv", "X", status);
    vvd(double(v7[1]), -0.1645663574562769, sp_tolerance, "sla::vxv", "Y", status);
    vvd(double(v7[2]), -0.5093390925544726, sp_tolerance, "sla::vxv", "Z", status);

    // do same tests in double precision
    vector<double> dav;
    matrix<double> drm1;
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

    matrix<double> drm2;
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

    matrix<double> drm;
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

    vector<double> dv1;
    dcs2c({3.0123, -0.999}, dv1);
    vvd(dv1[0], -0.5366267667260525, dp_tolerance, "sla::dcs2c", "X", status);
    vvd(dv1[1], 0.06977111097651444, dp_tolerance, "sla::dcs2c", "Y", status);
    vvd(dv1[2], -0.8409302618566215, dp_tolerance, "sla::dcs2c", "Z", status);

    vector<double> dv2, dv3;
    dmxv(drm1, dv1, dv2);
    dmxv(drm2, dv2, dv3);
    vvd(dv3[0], -0.7267487768696160, dp_tolerance, "sla::dmxv", "X", status);
    vvd(dv3[1], 0.5011537352639822, dp_tolerance, "sla::dmxv", "Y", status);
    vvd(dv3[2], 0.4697671220397141, dp_tolerance, "sla::dmxv", "Z", status);

    vector<double> dv4;
    dimxv(drm, dv3, dv4);
    vvd(dv4[0], -0.5366267667260526, dp_tolerance, "sla::dimxv", "X", status);
    vvd(dv4[1], 0.06977111097651445, dp_tolerance, "sla::dimxv", "Y", status);
    vvd(dv4[2], -0.8409302618566215, dp_tolerance, "sla::dimxv", "Z", status);

    vector<double> dv5;
    dm2av(drm, dv5);
    vvd(dv5[0], 0.006889040510209034, dp_tolerance, "sla::dm2av", "X", status);
    vvd(dv5[1], -1.577473205461961, dp_tolerance, "sla::dm2av", "Y", status);
    vvd(dv5[2], 0.5201843672856759, dp_tolerance, "sla::dm2av", "Z", status);

    for (int j = 0; j < 3; j++) {
        dv5[j] *= 1000.0;
    }
    vector<double> dv6;
    double dvm = dvn(dv5, dv6);
    vvd(dv6[0], 0.004147420704640065, dp_tolerance, "sla::dvn", "X", status);
    vvd(dv6[1], -0.9496888606842218, dp_tolerance, "sla::dvn", "Y", status);
    vvd(dv6[2], 0.3131674740355448, dp_tolerance, "sla::dvn", "Z", status);
    vvd(dvm, 1661.042127339937, 1.0e-9, "sla::DVN", "m", status);

    vvd(dvdv(dv6, dv1), -0.3318384698006295, dp_tolerance, "sla::DVN", " ", status);

    vector<double> dv7;
    dvxv(dv6, dv1, dv7);
    vvd(dv7[0], 0.7767720597123304, dp_tolerance, "sla::dvxv", "X", status);
    vvd(dv7[1], -0.1645663574562769, dp_tolerance, "sla::dvxv", "Y", status);
    vvd(dv7[2], -0.5093390925544726, dp_tolerance, "sla::dvxv", "Z", status);
}

/// Tests sla::bear(), sla::dbear(), sla::pav(), and sla::dpav() functions.
void t_bear(bool& status) {
    vector<float> fv1, fv2;
    vector<double> dv1, dv2;
    const double a1 = 1.234;
    const double b1 = -0.123;
    const double a2 = 2.345;
    const double b2 = 0.789;

    vvd(double(bear(float(a1), float(b1), float(a2), float(b2))), 0.7045970341781791, 1.0e-6, "sla::bear", " ", status);
    vvd(dbear(a1, b1, a2, b2), 0.7045970341781791, 1.0e-12, "sla::dbear", " ", status);
    dcs2c({a1, b1}, dv1);
    dcs2c({a2, b2}, dv2);

    for (int i = 0; i < 3; i++) {
        fv1[i] = float(dv1[i]);
        fv2[i] = float(dv2[i]);
    }
    vvd(double(pav(fv1, fv2 )), 0.7045970341781791, 1.0e-6, "sla::pav", " ", status);
    vvd(dpav(dv1, dv2), 0.7045970341781791, 1.0e-12, "sla::dpav", " ", status);
}

/// Tests sla::zd() function.
void t_zd(bool& status) {
    vvd(zd(-1.023, -0.876, -0.432), 0.8963914139430839, 1.0e-12, "sla::zd", " ", status);
}

///////////////////////////////////////////////////////////////////////////////
// MODULE ENTRY POINT
///////////////////////////////////////////////////////////////////////////////

/// Tests all SLALIB functions and procedures.
bool sla_test() {
    bool status = true;
    t_airmas(status);
    t_bear(status);
    t_vecmat(status);
    t_zd(status);
    return status;
}

}
