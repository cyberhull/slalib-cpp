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
 * Projects spherical coordinates onto tangent plane: "gnomonic" projection - "standard coordinates" (double
 * precision).
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param point Spherical coordinates of point to be projected.
 * @param tangent Spherical coordinates of tangent point.
 * @param xi Return value: rectangular coordinate on tangent plane.
 * @param eta Return value: rectangular coordinate on tangent plane.
 * @return A `TPPStatus` constant.
 */
TPPStatus ds2tp(const Spherical<double>& point, const Spherical<double>& tangent, double& xi, double& eta) {
    constexpr double TINY = 1e-6;

    // trig functions
    const double sin_tdec = std::sin(tangent.get_dec());
    const double sin_dec = std::sin(point.get_dec());
    const double cos_tdec = std::cos(tangent.get_dec());
    const double cos_dec = std::cos(point.get_dec());
    const double ra_diff = point.get_ra() - tangent.get_ra();
    const double sin_ra_diff = std::sin(ra_diff);
    const double cos_ra_diff = std::cos(ra_diff);

    // reciprocal of star vector length to tangent plane
    double denom = sin_dec * sin_tdec + cos_dec * cos_tdec * cos_ra_diff;

    // handle vectors too far from axis
    TPPStatus result;
    if (denom > TINY) {
        result = TPP_OK;
    } else if (denom >= 0.0) {
        result = TPP_TOO_FAR;
        denom = TINY;
    } else if (denom > -TINY) {
        result = TPP_ASTAR_ON_TP;
        denom = -TINY;
    } else {
        result = TPP_ASTAR_TOO_FAR;
    }

    // compute tangent plane coordinates (even in dubious cases)
    xi = cos_dec * sin_ra_diff / denom;
    eta = (sin_dec * cos_tdec - cos_dec * sin_tdec * cos_ra_diff) / denom;
    return result;
}

}
