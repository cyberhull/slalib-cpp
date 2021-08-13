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
    vector<float> fv1, fv2;
    vector<double> dv1, dv2;
    constexpr double a1 = 1.234;
    constexpr double b1 = -0.123;
    constexpr double a2 = 2.345;
    constexpr double b2 = 0.789;

    vvd(bear(float(a1), float(b1), float(a2), float(b2)), 0.7045970341781791, 1.0e-6, "sla::bear", " ", status);
    vvd(dbear(a1, b1, a2, b2), 0.7045970341781791, 1.0e-12, "sla::dbear", " ", status);
    dcs2c({a1, b1}, dv1);
    dcs2c({a2, b2}, dv2);

    for (int i = 0; i < 3; i++) {
        fv1[i] = float(dv1[i]);
        fv2[i] = float(dv2[i]);
    }
    vvd(pav(fv1, fv2 ), 0.7045970341781791, 1.0e-6, "sla::pav", " ", status);
    vvd(dpav(dv1, dv2), 0.7045970341781791, 1.0e-12, "sla::dpav", " ", status);
}

// tests sla::caf2r() and sla::daf2r() procedures
static void t_caf2r(bool& status) {
    float sp_radians;
    D2RStatus result = caf2r(76, 54, 32.1f, sp_radians);
    vvd(sp_radians, 1.342313819975276, 1.0e-6, "sla::caf2r", "R", status);
    viv(result, 0, "sla::caf2r", "S", status );

    double dp_radians;
    result = daf2r(76, 54, 32.1, dp_radians);
    vvd(dp_radians, 1.342313819975276, 1.0e-12, "sla::daf2r", "R", status);
    viv(result, 0, "sla::daf2r", "S", status);
}

// tests sla::caldj() procedure
static void t_caldj(bool& status) {
    double mjd;
    auto result = caldj(1999, 12, 31, mjd);
    vvd(mjd, 51543.0, 0.0, "sla::caldj", "D", status);
    viv(result, 0, "sla::caldj", "S", status);
}

// tests sla::clyd() and sla::calyd() procedures
static void t_calyd(bool& status) {
    int j_year, j_day;
    G2JStatus result;

    result = calyd(46, 4, 30, j_year, j_day);
    viv(j_year, 2046, "sla::calyd", "Y", status);
    viv(j_day, 120, "sla::calyd", "D", status);
    viv(result, 0, "sla::calyd", "S", status);
    result = clyd (-5000, 1, 1, j_year, j_day);
    viv(result, 1, "sla::clyd", "illegal year", status);
    result = clyd (1900, 0, 1, j_year, j_day);
    viv(result, 2, "sla::clyd", "illegal month", status);
    result = clyd (1900, 2, 29, j_year, j_day);
    viv(j_year, 1900, "sla::clyd", "illegal day (Y)", status);
    viv(j_day, 61, "sla::clyd", "illegal day (D)", status);
    viv(result, 3, "sla::clyd", "illegal day (S)", status);
    result = clyd (2000, 2, 29, j_year, j_day);
    viv(j_year, 2000, "sla::clyd", "Y", status);
    viv(j_day, 60, "sla::clyd", "D", status);
    viv(result, 0, "sla::clyd", "S", status);
}

// tests sla::cc2s() and dcc2s() procedures
static void t_cc2s(bool& status) {
    const vector<float> V = {100.0f, -50.0f, 25.0f};
    SphericalDir<float> sp_spherical;
    cc2s(V, sp_spherical);
    vvd(sp_spherical.sd_a, -0.4636476090008061, 1.0e-6, "sla::cc2s", "A", status );
    vvd(sp_spherical.sd_b, 0.2199879773954594, 1.0e-6, "sla::cc2s", "B", status );

    const vector<double> DV = {100.0, -50.0, 25.0};
    SphericalDir<double> dp_spherical;
    dcc2s(DV, dp_spherical);
    vvd(dp_spherical.sd_a, -0.4636476090008061, 1.0e-12, "sla::dcc2s", "A", status );
    vvd(dp_spherical.sd_b, 0.2199879773954594, 1.0e-12, "sla::dcc2s", "B", status );
}

// tests sla::cldj() procedure
static void t_cldj(bool& status) {
    double mjd;
    auto result = cldj(1899, 12, 31, mjd);
    vvd(mjd, 15019.0, 0.0, "sla::cldj", "D", status);
    viv(result, 0, "sla::cldj", "S", status);
}

// tests sla::e2h(), sla::de2h(), sla::h2e(), and sla::dh2e() procedures
static void t_e2h(bool& status) {
    double d_ha = -0.3;
    double d_dec = -1.1;
    double d_phi = -0.7;
    auto ha = float(d_ha);
    auto dec = float(d_dec);
    auto phi = float(d_phi);

    double d_azimuth, d_elevation;
    de2h(d_ha, d_dec, d_phi, d_azimuth, d_elevation);
    vvd(d_azimuth, 2.820087515852369, 1.0e-12, "sla::de2h", "Az", status);
    vvd(d_elevation, 1.132711866443304, 1.0e-12, "sla::de2h", "El", status);

    float azimuth, elevation;
    e2h(ha, dec, phi, azimuth, elevation );
    vvd(azimuth, 2.820087515852369, 1.0e-6, "sla::e2h", "Az", status);
    vvd(elevation, 1.132711866443304, 1.0e-6, "sla::e2h", "El", status);

    dh2e(d_azimuth, d_elevation, d_phi, d_ha, d_dec);
    vvd(d_ha, -0.3, 1.0e-12, "sla::dh2e", "HA", status);
    vvd(d_dec, -1.1, 1.0e-12, "sla::dh2e", "Dec", status);

    h2e(azimuth, elevation, phi, ha, dec);
    vvd(ha, -0.3, 1.0e-6, "sla::h2e", "HA", status);
    vvd(dec, -1.1, 1.0e-6, "sla::h2e", "Dec", status);
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
    vector<float> av;
    matrix<float> rm1;
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
    matrix<float> rm2;
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
    matrix<float> rm;
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
    vector<float> v1;
    cs2c({3.0123f, -0.999f}, v1 );
    vvd(v1[0], -0.5366267667260525, sp_tolerance, "sla::cs2c", "X", status);
    vvd(v1[1], 0.06977111097651444, sp_tolerance, "sla::cs2c", "Y", status);
    vvd(v1[2], -0.8409302618566215, sp_tolerance, "sla::cs2c", "Z", status);

    // rotate the vector using the two matrices sequentially
    vector<float> v2, v3;
    mxv(rm1, v1, v2);
    mxv(rm2, v2, v3);
    vvd(v3[0], -0.7267487768696160, sp_tolerance, "sla::mxv", "X", status);
    vvd(v3[1], 0.5011537352639822, sp_tolerance, "sla::mxv", "Y", status);
    vvd(v3[2], 0.4697671220397141, sp_tolerance, "sla::mxv", "Z", status);

    // de-rotate the vector using the combined matrix
    vector<float> v4;
    imxv(rm, v3, v4);
    vvd(v4[0], -0.5366267667260526, sp_tolerance, "sla::imxv", "X", status);
    vvd(v4[1], 0.06977111097651445, sp_tolerance, "sla::imxv", "Y", status);
    vvd(v4[2], -0.8409302618566215, sp_tolerance, "sla::imxv", "Z", status);

    // convert the combined matrix into an axial vector
    vector<float> v5;
    m2av(rm, v5);
    vvd(v5[0], 0.006889040510209034, sp_tolerance, "sla::m2av", "X", status);
    vvd(v5[1], -1.577473205461961, sp_tolerance, "sla::m2av", "Y", status);
    vvd(v5[2], 0.5201843672856759, sp_tolerance, "sla::m2av", "Z", status);

    // multiply the axial vector by a scalar and then normalize
    for (int i = 0; i < 3; i++) {
        v5[i] *= 1000.0f;
    }
    vector<float> v6;
    float vm = vn(v5, v6);
    vvd(v6[0], 0.004147420704640065, sp_tolerance, "sla::vn", "X", status);
    vvd(v6[1], -0.9496888606842218, sp_tolerance, "sla::vn", "Y", status);
    vvd(v6[2], 0.3131674740355448, sp_tolerance, "sla::vn", "Z", status);
    vvd(vm, 1661.042127339937, 1.0e-3, "sla::vn", "M", status);

    // calculate dot product with the original vector
    vvd(vdv(v6, v1), -0.3318384698006295, sp_tolerance, "sla::vn", " ", status);

    // calculate cross product with the original vector
    vector<float> v7;
    vxv(v6, v1, v7);
    vvd(v7[0], 0.7767720597123304, sp_tolerance, "sla::vxv", "X", status);
    vvd(v7[1], -0.1645663574562769, sp_tolerance, "sla::vxv", "Y", status);
    vvd(v7[2], -0.5093390925544726, sp_tolerance, "sla::vxv", "Z", status);

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
    vvd(dvm, 1661.042127339937, 1.0e-9, "sla::dvn", "M", status);

    vvd(dvdv(dv6, dv1), -0.3318384698006295, dp_tolerance, "sla::dvn", " ", status);

    vector<double> dv7;
    dvxv(dv6, dv1, dv7);
    vvd(dv7[0], 0.7767720597123304, dp_tolerance, "sla::dvxv", "X", status);
    vvd(dv7[1], -0.1645663574562769, dp_tolerance, "sla::dvxv", "Y", status);
    vvd(dv7[2], -0.5093390925544726, dp_tolerance, "sla::dvxv", "Z", status);
}

// tests sla::zd() function
static void t_zd(bool& status) {
    vvd(zd(-1.023, -0.876, -0.432), 0.8963914139430839, 1.0e-12, "sla::zd", " ", status);
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
    viv((int) result.get_sign(), int('-'), "sla::DR2TF", "S", status);
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

    double REFA2, REFB2;
    atmdsp(275.9, 709.3, 0.9, 0.77, refa, refb, 0.5, REFA2, REFB2);
    vvd(REFA2, 2.034523658888048e-4, 1.0e-12, "sla::atmdsp", "a", status);
    vvd(REFB2, -2.250855362179e-7, 1.0e-15, "sla::atmdsp", "b", status);

    const SphericalDir<double> spherical1 = {0.345, 0.456};
    vector<double> cartesian1, cartesian2;
    dcs2c(spherical1, cartesian1);
    refv(cartesian1, refa, refb, cartesian2);
    vvd(cartesian2[0], 0.8447487047790478, 1.0e-12, "sla::refv", "x1", status);
    vvd(cartesian2[1], 0.3035794890562339, 1.0e-12, "sla::refv", "y1", status);
    vvd(cartesian2[2], 0.4407256738589851, 1.0e-12, "sla::refv", "z1", status);

    const SphericalDir<double> spherical2 = {3.7, 0.03};
    dcs2c(spherical2, cartesian1);
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
    matrix<double> rm;
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
    t_cc2s(status);
    t_cldj(status);
    t_e2h(status);
    t_vecmat(status);
    t_zd(status);
    t_cd2tf(status);
    t_cr2af(status);
    t_cr2tf(status);
    t_ctf2d(status);
    t_ctf2r(status);
    t_range(status);
    t_ranorm(status);
    t_ref(status);
    t_ecmat(status);
    return status;
}

}
