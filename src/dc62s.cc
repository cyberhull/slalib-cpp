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
 * Converts and velocity in Cartesian coordinates to spherical coordinates (double precision).
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param cartesian Cartesian position & velocity vector.
 * @param spherical Longitude, latitude (radians), radial coordinate, longitude derivative, latitude derivative
 *   (radians per unit time), radial derivative.
 */
void dc62s(const VectorPV<double>& cartesian, SphericalPV<double>& spherical) {
    // components of the position/velocity vector
    double x = cartesian.get_x();
    double y = cartesian.get_y();
    double z = cartesian.get_z();
    const double xd = cartesian.get_dx();
    const double yd = cartesian.get_dy();
    const double zd = cartesian.get_dz();

    // component of dist in XY plane squared
    double rxy2 = x * x + y * y;

    // modulus squared
    double r2 = rxy2 + z * z;

    // protection against null vector
    if (r2 == 0.0) {
        x = xd;
        y = yd;
        z = zd;
        rxy2 = x * x + y * y;
        r2 = rxy2 + z * z;
    }

    // position and velocity in spherical coordinates
    const double rxy = std::sqrt(rxy2);
    const double xyp = x * xd + y * yd;
    if (rxy2 != 0.0) {
        spherical.set_longitude(std::atan2(y, x));
        spherical.set_latitude(std::atan2(z, rxy));
        spherical.set_dlong((x * yd - y * xd) / rxy2);
        spherical.set_dlat((zd * rxy2 - z * xyp) / (r2 * rxy));
    } else {
        spherical.set_longitude(0.0);
        spherical.set_latitude((z != 0.0)? std::atan2(z, rxy): 0.0);
        spherical.set_dlong(0.0);
        spherical.set_dlat(0.0);
    }
    const double dist = std::sqrt(r2);
    spherical.set_dist(dist);
    spherical.set_ddist((dist != 0.0)? (xyp + z * zd) / dist: 0.0);
}

}
