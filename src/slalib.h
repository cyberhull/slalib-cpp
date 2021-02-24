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
#ifndef _SLALIB_H
#define _SLALIB_H

namespace sla {

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

} // sla

#endif // _SLALIB_H
