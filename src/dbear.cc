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
 * Calculates bearing (position angle) of one point on a sphere relative to another given spherical coordinates of
 * the point (double precision).
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param da RA/longitude/etc. and Dec/latitude/etc. of the first point (radians).
 * @param db RA/longitude/etc. and Dec/latitude/etc. of the second point (radians).
 * @return Bearing (position angle), in radians (range +/- pi), of point `db` as seen from point `da`; if
 *  `db` is due east of `da` the bearing is +pi/2; if the two points are coincident, zero is returned.
 */
double dbear(const Spherical<double>& da, const Spherical<double>& db) {
    const double cos_b_dec = std::cos(db.get_dec());
    const double d_ra = db.get_ra() - da.get_ra();
    const double x = std::sin(db.get_dec()) * std::cos(da.get_dec()) - cos_b_dec * std::sin(da.get_dec()) * std::cos(d_ra);
    const double y = std::sin(d_ra) * cos_b_dec;
    return (x != 0.0 || y != 0.0)? std::atan2(y, x): 0.0;
}

}
