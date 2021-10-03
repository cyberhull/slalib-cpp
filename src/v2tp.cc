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
 * (single precision).
 *
 *  This function is the Cartesian equivalent of the function sla::s2tp().
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param point Direction cosines of star; if this vector is of zero length, the results will be wrong.
 * @param tangent Direction cosines of tangent point; if this vector is not of unit length, the results will be wrong; if
 *   this vector points at a pole, the returned `xi`,`eta` will be based on the arbitrary assumption that the RA of
 *   the tangent point is zero.
 * @param xi Return value: tangent plane coordinate of star.
 * @param eta Return value: tangent plane coordinate of star.
 * @return A `TPPStatus` constant (status).
 */
TPPStatus v2tp(const Vector<float> point, const Vector<float> tangent, float& xi, float& eta) {
    constexpr float TINY = 1E-6f;

    const float x = point[0];
    const float y = point[1];
    const float z = point[2];
    float x0 = tangent[0];
    const float y0 = tangent[1];
    const float z0 = tangent[2];
    const float r2 = x0 * x0 + y0 * y0;
    float r = std::sqrt(r2);
    if (r == 0.0f) {
        x0 = r = 1e-20f;
    }
    const float w = x * x0 + y * y0;
    float d = w + z * z0;
    TPPStatus status;
    if (d > TINY) {
        status = TPP_OK;
    } else if (d >= 0.0f) {
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
