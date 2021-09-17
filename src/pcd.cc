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

namespace sla {
/**
 * Applies pincushion/barrel distortion to a tangent-plane [x,y] (double precision).
 *
 * The distortion is of the form RP = R*(1 + C*R**2), where R is the radial distance from the tangent point, C is the
 * `disco` argument, and RP is the radial distance in the presence of the distortion.
 *
 * For pincushion distortion, C is positive; for barrel distortion, C is negative.
 *
 * For `x`, `y` in units of one projection radius (in the case of a camera matrix, the focal length), the following
 * `disco` values apply:
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
 * There is a companion function, sla::unpcd(), which performs the inverse operation.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param disco Pincushion/barrel distortion coefficient.
 * @param x Input/return value: tangent-plane X coordinate (input) / distorted X coordinate (returned).
 * @param y Input/return value: tangent-plane Y coordinate (input) / distorted Y coordinate (returned).
 */
void pcd(double disco, double& x, double& y) {
    const double f = 1.0 + disco * (x * x + y * y);
    x *= f;
    y *= f;
}

}
