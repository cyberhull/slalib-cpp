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
 * Given the direction cosines of a star and of the tangent point, determines the star's tangent-plane coordinates
 * (double precision).
 *
 *  This function is the Cartesian equivalent of the function sla::ds2tp().
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param v Direction cosines of star; if this vector is of zero length, the results will be wrong.
 * @param v0 Direction cosines of tangent point; if this vector is not of unit length, the results will be wrong; if
 *   this vector points at a pole, the returned `xi`,`eta` will be based on the arbitrary assumption that the RA of
 *   the tangent point is zero.
 * @param xi Return value: tangent plane coordinate of star.
 * @param eta Return value: tangent plane coordinate of star.
 * @return A `TPPStatus` constant (status).
 */
TPPStatus dv2tp(const Vector<double> v, const Vector<double> v0, double& xi, double& eta) {
    constexpr double TINY = 1e-6;

    const double x = v[0];
    const double y = v[1];
    const double z = v[2];
    double x0 = v0[0];
    const double y0 = v0[1];
    const double z0 = v0[2];
    const double r2 = x0 * x0 + y0 * y0;
    double r = std::sqrt(r2);
    if (r == 0.0f) {
        x0 = r = 1e-20;
    }
    const double w = x * x0 + y * y0;
    double d = w + z * z0;
    TPPStatus status;
    if (d > TINY) {
        status = TPP_OK;
    } else if (d >= 0.0) {
        status = TPP_TOO_FAR;
        d = TINY;
    } else if (d > -TINY) {
        status = TPP_ASTAR_ON_TP;
        d = -TINY;
    } else {
        status = TPP_ASTAR_TOO_FAR;
    }
    d *= r;
    xi = (y * x0 - x * y0) / d;
    eta = (z * r2 - z0 * w) / d;
    return status;
}

}
