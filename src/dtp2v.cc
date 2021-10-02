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
 * Given the tangent-plane coordinates of a star and the direction cosines of the tangent point, determines the
 * direction cosines of the star (double precision).
 *
 * This function is the Cartesian equivalent of the function sla::tp2s().
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param xi Tangent plane coordinate of the star.
 * @param eta Tangent plane coordinate of the star.
 * @param v0 Direction cosines of the tangent point; if this vector is not of unit length, the returned vector `v`
 *   will be wrong; if this vector points at a pole, the returned vector `v` will be based on the arbitrary
 *   assumption that the RA of the tangent point is zero.
 * @param v Return value: direction cosines of star.
 */
void dtp2v(double xi, double eta, const Vector<double> v0, Vector<double> v) {
    double x = v0[0];
    const double y = v0[1];
    const double z = v0[2];
    const double f = std::sqrt(1.0 + xi * xi + eta * eta);
    double r = std::sqrt(x * x + y * y);
    if (r == 0.0) {
        x = r = 1e-20;
    }
    v[0] = (x - (xi * y + eta * x * z) / r) / f;
    v[1] = (y + (xi * x - eta * y * z) / r) / f;
    v[2] = (z + eta * r) / f;
}

}
