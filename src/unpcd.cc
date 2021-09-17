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
#include "f77_utils.h"

namespace sla {

/**
 * Removes pincushion/barrel distortion from a distorted [x,y] to give tangent-plane [x,y] (double precision).
 *
 * The distortion is of the form RP = R*(1+C*R^2), where R is the radial distance from the tangent point, C is the
 * `disco` argument, and RP is the radial distance in the presence of the distortion.
 *
 * For pincushion distortion, C is positive; for barrel distortion, C is negative.
 *
 * For `x`, `y` in "radians" - units of one projection radius, which in the case of a photograph is the focal
 * length of the camera - the following `disco` values apply:
 *
 *   Geometry        `disco`
 *   --------------  ---------
 *   astrograph         0.0
 *   Schmidt           -0.3333
 *   AAT PF doublet  +147.069
 *   AAT PF triplet  +178.585
 *   AAT f/8          +21.20
 *   JKT f/8          +13.32
 *
 * The present function is a rigorous inverse of the companion function sla::pcd(). The expression for RP above is
 * rewritten in the form x^3+a*x+b=0 and solved by standard techniques.
 *
 * Cases where the cubic has multiple real roots can sometimes occur, corresponding to extreme instances of barrel
 * distortion where up to three different undistorted [X,Y]s all produce the same distorted [x,y]. However, only one
 * solution is returned, the one that produces the smallest change in [x,y].
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param disco Pincushion/barrel distortion coefficient.
 * @param x Input/return value: distorted X coordinate (input) / tangent-plane X coordinate (returned).
 * @param y Input/return value: distorted Y coordinate (input) / tangent-plane Y coordinate (returned).
 */
void unpcd(double disco, double& x, double& y) {
    constexpr double ONE_THIRD = 1.0 / 3.0;
    constexpr double PI2 = 6.283185307179586476925286766559;

    // distance of the point from the origin
    const double rp = std::sqrt(x * x + y * y);

    // if zero, or if no distortion, no action is necessary
    if (rp != 0.0 && disco != 0.0) {

        // begin algebraic solution
        const double q = 1.0 / (3.0 * disco);
        const double r = rp / (2.0 * disco);
        double w = q * q * q + r * r;

        // continue if one real root, or three of which only one is positive
        double f;
        if (w >= 0.0) {
            const double d = std::sqrt(w);
            w = r + d;
            const double s1 = f_sign(std::pow(std::abs(w), ONE_THIRD), w);
            w = r - d;
            const double t = f_sign(std::pow(std::abs(w), ONE_THIRD), w);
            f = s1 + t;
        } else {
            // three different real roots: use geometrical method instead
            w = 2.0 / std::sqrt(-3.0 * disco);
            const double c = 4.0 * rp / (disco * w * w * w);
            const double s2 = std::sqrt(1.0 - std::min(c * c, 1.0));
            const double t3 = std::atan2(s2, c);

            // the three solutions
            const double f1 = w * std::cos((PI2 - t3) / 3.0);
            const double f2 = w * std::cos((t3) / 3.0);
            const double f3 = w * std::cos((PI2 + t3) / 3.0);

            // pick the one that moves [x,y] least
            const double w1 = std::abs(f1 - rp);
            const double w2 = std::abs(f2 - rp);
            const double w3 = std::abs(f3 - rp);
            if (w1 < w2) {
                if (w1 < w3) {
                    f = f1;
                } else {
                    f = f3;
                }
            } else {
                if (w2 < w3) {
                    f = f2;
                } else {
                    f = f3;
                }
            }

        }
        // remove the distortion
        f = f / rp;
        x = f * x;
        y = f * y;
    }
}

}
