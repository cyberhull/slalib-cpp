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
#ifndef SLALIB_H_INCLUDED
#define SLALIB_H_INCLUDED

namespace sla {

/// Status codes for the caf2r() procedure (degrees to radians conversion).
enum D2RStatus {
    D2R_OK = 0,         ///< all arguments fit their ranges, conversion successful
    D2R_BAD_DEGREES,    ///< degrees outside of range [0..359]
    D2R_BAD_ARCMINUTES, ///< minutes outside of range [0..59]
    D2R_BAD_ARCSECONDS  ///< seconds outside of range [0..60)
};

/// Status codes for the cldj(), caldj(), clyd(), and calyd() procedures (Gregorian to Julian calendar conversions).
enum G2JStatus {
    G2J_OK = 0,    ///< all arguments fit their ranges, conversion successful
    G2J_BAD_YEAR,  ///< year earlier (less) than -4699, output value(s) not calculated/returned
    G2J_BAD_MONTH, ///< month outside of range [1..12], output value(s) not calculated/returned
    G2J_BAD_DAY    ///< day outside [0..<days-in-given-month>] range, BUT output value(s) were calculated and returned
};

/// Generic 3-component vector of floating-point elements.
template<typename T>
using vector = T[3];

/// Generic 3x3 matrix of floating-point elements.
template<typename T>
using matrix = T[3][3];

/// Structure representing partial spherical coordinates: longitude/latitude, or right ascension/declination
template <typename T>
struct SphericalDir {
    T sd_a; ///< longitude or RA (radians)
    T sd_b; ///< latitude or Dec (radians)
};

/// Structure representing full spherical coordinates: longitude, latitude, and distance
template <typename T>
struct SphericalCoords {
    T sc_long; ///< longitude (radians)
    T sc_lat;  ///< latitude (radians)
    T sc_dist; ///< distance along long/lat ray
};

/**
 * Class representing various conversion results: days to hours, minutes, seconds; or radians to degrees, arcminutes,
 * arcseconds; etc. The same data structure has to be passed between routines interpreting it quite differently, hence
 * this implementation, having different interfaces to the same underlying data structure.
 */
class ConversionResult {
    int  cr_data[4]; ///< hours/minutes/seconds/fraction or degrees/arcminutes/arcseconds/fraction
    bool cr_sign;    ///< `false` for '-', `true` for '+'

public:
    // interface for radians/days to hours, minutes, seconds conversions
    [[nodiscard]] int get_hours() const { return cr_data[0]; }
    void set_hours(int hours) { cr_data[0] = hours; }
    [[nodiscard]] int get_minutes() const { return cr_data[1]; }
    void set_minutes(int minutes) { cr_data[1] = minutes; }
    [[nodiscard]] int get_seconds() const { return cr_data[2]; }
    void set_seconds(int seconds) { cr_data[2] = seconds; }

    // interface for radians to degrees, arcminutes, arcseconds conversion
    [[nodiscard]] int get_degrees() const { return cr_data[0]; }
    void set_degrees(int degrees) { cr_data[0] = degrees; }
    [[nodiscard]] int get_arcminutes() const { return cr_data[1]; }
    void set_arcminutes(int arcminutes) { cr_data[1] = arcminutes; }
    [[nodiscard]] int get_arcseconds() const { return cr_data[2]; }
    void set_arcseconds(int arcseconds) { cr_data[2] = arcseconds; }

    // interface shared by all conversions
    [[nodiscard]] int get_fraction() const { return cr_data[3]; }
    void set_fraction(int fraction) { cr_data[3] = fraction; }
    [[nodiscard]] char get_sign() const { return cr_sign? '+': '-'; }
    void set_sign(char sign) { cr_sign = sign == '+'; }
    void set_sign(bool sign) { cr_sign = sign; }
};

// auxiliary functions
int process_year_defaults(int year);
G2JStatus validate_gregorian_day(int year, int month, int day);

// library API
double airmas(double zenith_dist);
void av2m(const vector<float> vec, matrix<float> mat);
void dav2m(const vector<double> vec, matrix<double> mat);
void cc2s(const vector<float> cartesian, SphericalDir<float>& spherical);
void dcc2s(const vector<double> cartesian, SphericalDir<double>& spherical);
void cs2c(const SphericalDir<float>& spherical, vector<float> cartesian);
void dcs2c(const SphericalDir<double>& spherical, vector<double> cartesian);
void euler(const char* order, float phi, float theta, float psi, matrix<float> rmat);
void deuler(const char* order, double phi, double theta, double psi, matrix<double> rmat);
void imxv(const matrix<float> rm, const vector<float> va, vector<float> vb);
void dimxv(const matrix<double> rm, const vector<double> va, vector<double> vb);
void m2av(const matrix<float> rmat, vector<float> axis);
void dm2av(const matrix<double> rmat, vector<double> axis);
void mxm(const matrix<float> a, const matrix<float> b, matrix<float> c);
void dmxm(const matrix<double> a, const matrix<double> b, matrix<double> c);
void mxv(const matrix<float> rm, const vector<float> va, vector<float> vb);
void dmxv(const matrix<double> rm, const vector<double> va, vector<double> vb);
float vdv(const vector<float> va, const vector<float> vb);
double dvdv(const vector<double> va, const vector<double> vb);
float vn(const vector<float> v, vector<float> uv);
double dvn(const vector<double> v, vector<double> uv);
void vxv(const vector<float> va, const vector<float> vb, vector<float> vc);
void dvxv(const vector<double> va, const vector<double> vb, vector<double> vc);
double zd(double ha, double dec, double phi);
double bear(float a1, float b1, float a2, float b2);
double dbear(double a1, double b1, double a2, double b2);
float pav(const vector<float> v1, const vector<float> v2);
double dpav(const vector<double> v1, const vector<double> v2);
void e2h(float ha, float dec, float phi, float& azimuth, float& elevation);
void de2h(double ha, double dec, double phi, double& azimuth, double& elevation);
void h2e(float azimuth, float elevation, float phi, float& ha, float& dec);
void dh2e(double azimuth, double elevation, double phi, double& ha, double& dec);
D2RStatus caf2r(int degrees, int minutes, float seconds, float& radians);
D2RStatus daf2r(int degrees, int minutes, double seconds, double& radians);
G2JStatus cldj(int year, int month, int day, double& mjd);
G2JStatus caldj(int year, int month, int day, double& mjd);
G2JStatus clyd(int year, int month, int day, int& jyear, int& jday);
G2JStatus calyd(int year, int month, int day, int& jyear, int& jday);
void cd2tf(int ndp, float days, ConversionResult& result);
void dd2tf(int ndp, double days, ConversionResult& result);

} // sla

#endif // SLALIB_H_INCLUDED
