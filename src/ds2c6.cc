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
 * Converts position and velocity in spherical coordinates to Cartesian coordinates (double precision).
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param spv Position and velocity in spherical coordinates (radians and radians per unit time, respectively).
 * @param pv Cartesian position and velocity vector.
 */
void ds2c6(const SphericalPV<double>& spv, VectorPV<double> pv) {
    // useful functions
    const double sin_long = std::sin(spv.get_longitude());
    const double cos_long = std::cos(spv.get_longitude());
    const double sin_lat = std::sin(spv.get_latitude());
    const double cos_lat = std::cos(spv.get_latitude());
    const double dist_cos_lat = spv.get_dist() * cos_lat;
    const double x = dist_cos_lat * cos_long;
    const double y = dist_cos_lat * sin_long;
    const double dist_dlat = spv.get_dist() * spv.get_dlat();
    const double w = dist_dlat * sin_lat - cos_lat * spv.get_ddist();

    // position
    pv.set_x(x);
    pv.set_y(y);
    pv.set_z(spv.get_dist() * sin_lat);

    // velocity
    pv.set_dx(-y * spv.get_dlong() - w * cos_long);
    pv.set_dy(x * spv.get_dlong() - w * sin_long);
    pv.set_dz(dist_dlat * cos_lat + sin_lat * spv.get_ddist());
}

}
