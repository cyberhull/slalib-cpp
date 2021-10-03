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
 * Given the tangent-plane coordinates of a star and its direction cosines, determines the direction cosines of the
 * tangent-point (single precision).
 *
 * Cases where there is no solution can only arise near the poles. For example, it is clearly impossible for a star
 * at the pole itself to have a non-zero `xi` value, and hence it is meaningless to ask where the tangent point would
 * have to be.
 *
 * Also near the poles, cases can arise where there are two useful solutions. The return value indicates whether the
 * second of the two solutions returned is useful. Returned value of 1 indicates only one useful solution, the usual
 * case; under these circumstances, the second solution can be regarded as valid if the vector `solution2` is interpreted
 * as the "over-the-pole" case.
 *
 * This function is the Cartesian equivalent of the function sla::tps2c().
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param xi Tangent plane coordinate of star.
 * @param eta Tangent plane coordinate of star.
 * @param point Direction cosines of star; must be of unit length or the result will be wrong.
 * @param solution1 Return value: direction cosines of tangent point, solution 1.
 * @param solution2 Return value: direction cosines of tangent point, solution 2.
 * @return Number of solutions; 0: no solutions returned, 1: only the first solution is useful, 2: both solutions
 *   are useful.
 */
int tpv2c(float xi, float eta, const Vector<float> point, Vector<float> solution1, Vector<float> solution2) {
    const float x = point[0];
    const float y = point[1];
    const float z = point[2];
    const float rxy2 = x * x + y * y;
    const float xi2 = xi * xi;
    const float eta2p1 = eta * eta + 1.0f;
    const float sdf = z * std::sqrt(xi2 + eta2p1);
    const float r2 = rxy2 * eta2p1 - z * z * xi2;
    if (r2 > 0.0f) {
        float r = std::sqrt(r2);
        float c = (sdf * eta + r) / (eta2p1 * std::sqrt(rxy2 * (r2 + xi2)));
        solution1[0] = c * (x * r + y * xi);
        solution1[1] = c * (y * r - x * xi);
        solution1[2] = (sdf - eta * r) / eta2p1;
        r = -r;
        c = (sdf * eta + r) / (eta2p1 * std::sqrt(rxy2 * (r2 + xi2)));
        solution2[0] = c * (x * r + y * xi);
        solution2[1] = c * (y * r - x * xi);
        solution2[2] = (sdf - eta * r) / eta2p1;
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
