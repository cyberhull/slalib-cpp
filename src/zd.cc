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
 * Calculates zenith distance given HA, Dec, and observatory latitude (double precision).
 *
 * In some applications it will be important to specify the correct type of hour angle and declination in order to
 * produce the required type of zenith distance. In particular, it may be important to distinguish between the
 * zenith distance as affected by refraction, which would require the "observed" HA,Dec, and the zenith distance
 * in vacuo, which would require the "topocentric" HA, Dec. If the effects of diurnal aberration can be neglected,
 * the "apparent" HA, Dec may be used instead of the topocentric HA, Dec.
 *
 * In applications which involve many zenith distance calculations, rather than calling this function it will be more
 * efficient to use inline code, having previously computed fixed terms such as sine and cosine of the latitude, and
 * perhaps sine and cosine of the declination.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param dir Hour Angle and Declination (radians); not range-checked.
 * @param phi Observatory latitude (radians); not range-checked; must be geodetic (i.e. be the angle between the
 *   equator and the normal to the ellipsoid approximating the true Earth as defined by the Geodetic Reference System
 *   1980 (GRS-80)); in critical applications, corrections for polar motion should be applied.
 * @return Zenith distance, 0 to Pi.
 */
 // double ha, double dec
double zd(const Spherical<double>& dir, double phi) {
    const double sin_ha = std::sin(dir.get_ha());
    const double cos_ha = std::cos(dir.get_ha());
    const double sin_dec = std::sin(dir.get_dec());
    const double cos_dec = std::cos(dir.get_dec());
    const double sin_phi = std::sin(phi);
    const double cos_phi = std::cos(phi);
    const double x = cos_ha * cos_dec * sin_phi - sin_dec * cos_phi;
    const double y = sin_ha * cos_dec;
    const double z = cos_ha * cos_dec * cos_phi + sin_dec * sin_phi;
    return std::atan2(std::sqrt(x * x + y * y), z);
}

}
