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
 * Determines the RA,Dec of the tangent point from the tangent plane coordinates of a star of known RA,Dec (single
 * precision).
 *
 * The RA fields of the solutions are returned in the range 0-2pi; the Dec values of the solutions are returned in the
 * range +/-pi, but in the usual, non-pole-crossing, case, the range is +/-pi/2.
 *
 * Cases where there is no solution can only arise near the poles. For example, it is clearly impossible for a star
 * at the pole itself to have a non-zero `xi` value, and hence it is meaningless to ask where the tangent point would
 * have to be to bring about this combination of `xi` and Dec.
 *
 * Also near the poles, cases can arise where there are two useful solutions. The return value indicates whether the
 * second of the two solutions returned is useful. The return value of 1 indicates only one useful solution, the
 * usual case; under these circumstances, the second solution corresponds to the "over-the-pole" case, and this is
 * reflected in the RA and Dec values of `solution2` which are returned.
 *
 * This function is the spherical equivalent of the function sla::dtpv2c().
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param xi Tangent plane rectangular coordinate.
 * @param eta Tangent plane rectangular coordinate.
 * @param point Spherical coordinates of a star.
 * @param solution1 Return value: spherical coordinates of tangent point, solution 1.
 * @param solution2 Return value: spherical coordinates of tangent point, solution 2.
 * @return Number of solutions: 0 (no solutions returned), 1 (only the first solution is useful), or 2 (both
 *   solutions are useful).
 */
int tps2c(float xi, float eta, const Spherical<float>& point,
    Spherical<float>& solution1, Spherical<float>& solution2) {
    const float x2 = xi * xi;
    const float y2 = eta * eta;
    const float sin_dec = std::sin(point.get_dec());
    const float cos_dec = std::cos(point.get_dec());
    const float sdf = sin_dec * std::sqrt(1.0f + x2 + y2);
    const float r2 = cos_dec * cos_dec * (1.0f + y2) - sin_dec * sin_dec * x2;
    if (r2 >= 0.0f) {
        float r = std::sqrt(r2);
        float s = sdf - eta * r;
        float c = sdf * eta + r;
        if (xi == 0.0f && r == 0.0f) {
            r = 1.0f;
        }
        solution1.set_ra(ranorm(point.get_ra() - std::atan2(xi, r)));
        solution1.set_dec(std::atan2(s, c));
        r = -r;
        s = sdf - eta * r;
        c = sdf * eta + r;
        solution2.set_ra(ranorm(point.get_ra() - std::atan2(xi, r)));
        solution2.set_dec(std::atan2(s, c));
        if (std::abs(sdf) < 1.0f) {
            return 1;
        } else {
            return 2;
        }
    } else {
        return 0;
    }
}

}
