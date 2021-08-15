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
 * Calculates positions, velocities and accelerations for an altazimuth telescope mount (double precision).
 *
 * The velocities and accelerations assume constant declination and constant rate of change of hour angle (as for
 * tracking a star); to convert return values into practical degree- and second-based units:
 *
 *   angles * 360/2pi -> degrees
 *   velocities * (2pi/86400)*(360/2pi) -> degree/sec
 *   accelerations * ((2pi/86400)**2)*(360/2pi) -> degree/sec/sec
 *
 * Note that the seconds here are sidereal rather than SI. One sidereal second is about 0.99727 SI seconds.
 *
 * The velocity and acceleration factors assume the sidereal tracking case. Their respective numerical values
 * are (exactly) 1/240 and (approximately) 1/3300236.9.
 *
 * Refraction and deficiencies in the telescope mounting are ignored. The purpose of the routine is to give the
 * general form of the quantities. The details of a real telescope could profoundly change the results, especially
 * close to the zenith.
 *
 * In applications which involve many such calculations, rather than calling the present routine it will be more
 * efficient to use inline code, having previously computed fixed terms such as sine and cosine of latitude, and (for
 * tracking a star) sine and cosine of declination.
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param ha Hour angle (radians); topocentric; not range-checked.
 * @param dec Declination (radians); topocentric; not range-checked.
 * @param phi Observatory latitude (radians); geodetic as opposed to geocentric; not range-checked.
 * @param az Return value: azimuth (radians, range: [0..2pi]); north is zero, and east is +pi/2.
 * @param az_vel Return value: azimuth velocity (radians per radian of `ha`).
 * @param az_acc Return value: azimuth acceleration (radians per radian of `ha` squared).
 * @param el Return value: elevation (radians, range: [-pi..pi]).
 * @param el_vel Return value: elevation velocity (radians per radian of `ha`).
 * @param el_acc Return value: elevation acceleration (radians per radian of `ha` squared).
 * @param pa Return value: parallactic angle (radians, range: [-pi..pi]); +ve for a star west of the meridian and
 *   is the angle NP-star-zenith.
 * @param pa_vel Return value: parallactic angle velocity (radians per radian of `ha`).
 * @param pa_acc Return value: parallactic angle acceleration (radians per radian of `ha` squared).
 */
void altaz(double ha, double dec, double phi,
    double& az, double& az_vel, double& az_acc,
    double& el, double& el_vel, double& el_acc,
    double& pa, double& pa_vel, double& pa_acc) {

    constexpr double PI = 3.1415926535897932384626433832795;
    constexpr double PI2 = 6.283185307179586476925286766559;
    constexpr double EPSILON = 1.0e-30;

    // useful functions
    const double sin_ha = std::sin(ha);
    const double cos_ha = std::cos(ha);
    const double sin_dec = std::sin(dec);
    const double cos_dec = std::cos(dec);
    const double sin_phi = std::sin(phi);
    const double cos_phi = std::cos(phi);
    const double ch_cd = cos_ha * cos_dec;
    const double sd_cp = sin_dec * cos_phi;
    const double x = -ch_cd * sin_phi + sd_cp;
    const double y = -sin_ha * cos_dec;
    const double z = ch_cd * cos_phi + sin_dec * sin_phi;
    double r_squared = x * x + y * y;
    double r = std::sqrt(r_squared);

    // azimuth and elevation
    double azimuth = r_squared == 0.0? 0.0: std::atan2(y, x);
    if (azimuth < 0.0) {
        azimuth += PI2;
    }
    const double elevation = std::atan2(z, r);

    // parallactic angle
    const double c = cos_dec * sin_phi - cos_ha * sd_cp;
    const double s = sin_ha * cos_phi;
    const double p_angle = c * c + s * s > 0.0? std::atan2(s, c): PI - ha;

    // velocities and accelerations (clamped at zenith/nadir)
    if (r_squared < EPSILON) {
        r_squared = EPSILON;
        r = std::sqrt(r_squared);
    }
    const double p_vel = -x * cos_phi / r_squared;
    const double a_vel = sin_phi + z * p_vel;
    const double e_vel = cos_phi * y / r;
    const double edr = e_vel / r;
    const double a_acc = edr * (z * sin_phi + (2.0 - r_squared) * p_vel);
    const double e_acc = -r * p_vel * a_vel;
    const double p_acc = edr * (sin_phi + 2.0 * z * p_vel);

    // return results
    az = azimuth;
    az_vel = a_vel;
    az_acc = a_acc;
    el = elevation;
    el_vel = e_vel;
    el_acc = e_acc;
    pa = p_angle;
    pa_vel = p_vel;
    pa_acc = p_acc;
}

}
