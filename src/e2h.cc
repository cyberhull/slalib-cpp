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
 * Converts equatorial coordinates to horizon coordinates: HA,Dec to Az,El (single precision).
 *
 * In some applications, it will be important to specify the correct type of hour angle and declination in order to
 * produce the required type of azimuth and elevation. In particular, it may be important to distinguish between
 * elevation as affected by refraction, which would require the "observed" HA,Dec, and the elevation in vacuo, which
 * would require the "topocentric" (having origin at the observer) HA,Dec. If the effects of diurnal aberration can
 * be neglected, the "apparent" HA,Dec may be used instead of the topocentric HA,Dec.
 *
 * In applications which involve many such calculations, rather than calling this function it will be more efficient
 * to use inline code, having previously computed fixed terms such as sine and cosine of latitude, and (for tracking
 * a star) sine and cosine of declination.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param dir Hour angle and Declination (radians); not range-checked.
 * @param phi Observatory latitude, in radians; not range-checked; must be geodetic (i.e. be the angle between the
 *   equator and the normal to the ellipsoid approximating the true Earth as defined by the Geodetic Reference System
 *   1980 (GRS-80))); in critical applications, corrections for polar motion should be applied.
 * @param azimuth Output: azimuth; returned in the range 0-2Pi; north is zero, and east is +Pi/2.
 * @param elevation Output: elevation; returned in the range +/-Pi/2.
 */
void e2h(const Spherical<float>& dir, float phi, float& azimuth, float& elevation) {
    // useful trig functions
    const float sin_ha = std::sin(dir.get_ha());
    const float cos_ha = std::cos(dir.get_ha());
    const float sin_dec = std::sin(dir.get_dec());
    const float cos_dec = std::cos(dir.get_dec());
    const float sin_phi = std::sin(phi);
    const float cos_phi = std::cos(phi);

    // Az,El as x,y,z
    const float x = -cos_ha * cos_dec * sin_phi + sin_dec * cos_phi;
    const float y = -sin_ha * cos_dec;
    const float z = cos_ha * cos_dec * cos_phi + sin_dec * sin_phi;

    // to spherical coordinates
    float a;
    const float r = std::sqrt(x * x + y * y);
    if (r == 0.0f) {
        a = 0.0f;
    } else {
        a = std::atan2(y, x);
        if (a < 0.0f) {
            constexpr float R2PI = 6.283185307179586476925286766559f;
            a += R2PI;
        }
    }
    azimuth = a;
    elevation = std::atan2(z, r);
}

}
