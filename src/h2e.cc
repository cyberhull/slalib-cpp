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
 * Converts horizon coordinates to equatorial coordinates: Az,El to HA,Dec (single precision).
 *
 * In some applications it will be important to specify the correct type of elevation in order to produce the
 * required type of HA,Dec. In particular, it may be important to distinguish between the elevation as affected by
 * refraction, which will yield the "observed" HA,Dec, and the elevation in vacuo, which will yield the "topocentric"
 * HA,Dec. If the effects of diurnal aberration can be neglected, the topocentric HA,Dec may be used as an
 * approximation to the "apparent" HA,Dec.
 *
 * In applications which involve many such calculations, rather than calling the present routine it will be more
 * efficient to use inline code, having previously computed fixed terms such as sine and cosine of the latitude.
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param azimuth Azimuth, in radians; not range-checked; the sign convention for azimuth is north zero, east +Pi/2.
 * @param elevation Elevation, in radians; not range-checked.
 * @param phi Observatory latitude, in radians; not range-checked; the latitude is (in principle) geodetic; in
 *   critical applications, corrections for polar motion should be applied.
 * @param ha Return value: hour angle (range: +/-pi) and declination (range: +/-pi/2).
 */
void h2e(float azimuth, float elevation, float phi, Spherical<float>& dir) {
    // useful trig functions.
    const double sin_az = std::sin(azimuth);
    const double cos_az = std::cos(azimuth);
    const double sin_el = std::sin(elevation);
    const double cos_el = std::cos(elevation);
    const double sin_phi = std::sin(phi);
    const double cos_phi = std::cos(phi);

    // Az,El,Phi as x,y,z.
    const double x = -cos_az * cos_el * sin_phi + sin_el * cos_phi;
    const double y = -sin_az * cos_el;
    const double z = cos_az * cos_el * cos_phi + sin_el * sin_phi;

    // to ha,Dec
    const double r = std::sqrt(x * x + y * y);
    dir.set_ha((r == 0.0)? 0.0f: float(std::atan2(y, x)));
    dir.set_dec((float) std::atan2(z, r));
}

}
