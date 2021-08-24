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
 * Converts HA, Dec to Parallactic Angle (double precision).
 *
 * The parallactic angle at a point in the sky is the position angle of the vertical, i.e. the angle between the
 * direction to the pole and to the zenith. In precise applications care must be taken only to use geocentric apparent
 * HA,Dec and to consider separately the effects of atmospheric refraction and telescope mount errors.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param dir Hour angle and declination (radians; geocentric apparent).
 * @param phi Observatory latitude (radians;  geodetic).
 * @return Parallactic angle; at the pole, zero is returned.
 */
double pa(const Spherical<double>& dir, double phi) {
    const double cp = std::cos(phi);
    const double sqsz = cp * std::sin(dir.get_ha());
    double cqsz = std::sin(phi) * std::cos(dir.get_dec()) - cp * std::sin(dir.get_dec()) * std::cos(dir.get_ha());
    if (sqsz == 0.0 && cqsz == 0.0) {
        cqsz = 1.0;
    }
    return std::atan2(sqsz, cqsz);
}

}
