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
 * Adjusts an unrefracted Cartesian vector to include the effect of atmospheric refraction, using the simple
 * A tan Z + B tan**3 Z model.
 *
 * This function applies the adjustment for refraction in the opposite sense to the usual one -- it takes an
 * unrefracted (in vacuo) position and produces an observed (refracted) position, whereas the A tan Z + B tan**3 Z
 * model strictly applies to the case where an observed position is to have the refraction removed. The unrefracted
 * to refracted case is harder, and requires an inverted form of the text-book refraction models; the algorithm used
 * here is equivalent to one iteration of the Newton-Raphson method applied to the above formula.
 *
 * Though optimized for speed rather than precision, the present routine achieves consistency with the refracted-
 * to-unrefracted A tan Z + B tan**3 Z model at better than 1 microarcsecond within 30 degrees of the zenith and
 * remains within 1 milliarcsecond to beyond ZD 70 degrees. The inherent accuracy of the model is, of course, far
 * worse than this - see the documentation for sla::refco() for more information.
 *
 * At low elevations (below about 3 degrees), the refraction correction is held back to prevent arithmetic problems
 * and wildly wrong results. For optical/IR wavelengths, over a wide range of observer heights and corresponding
 * temperatures and pressures, the following levels of accuracy (arcseconds, worst case) are achieved, relative to
 * numerical integration through a model atmosphere:
 *
 *   ZD    error
 *   80      0.7
 *   81      1.3
 *   82      2.5
 *   83      5
 *   84     10
 *   85     20
 *   86     55
 *   87    160
 *   88    360
 *   89    640
 *   90   1100
 *   91   1700   } relevant only to
 *   92   2600   } high-elevation sites
 *
 * The results for radio are slightly worse over most of the range, becoming significantly worse below ZD==88 and
 * unusable beyond ZD==90.
 *
 * See also the function sla::refz(), which performs the adjustment to the zenith distance rather than in Cartesian
 * Az/El coordinates. The present function is faster than sla::refz() and, except very low down, is equally accurate
 * for all practical purposes. However, beyond about ZD==84 degrees sla::refz() should be used, and for the utmost
 * accuracy iterative use of sla::refro() should be considered.
 *
 * Original FORTRAN code by Rutherford Appleton Laboratory / P.T. Wallace.
 *
 * @param vu Unrefracted position of the source (Az/El 3-component vector).
 * @param refa Tan Z coefficient (radians).
 * @param refb Tan**3 Z coefficient (radians).
 * @param vr Return value: refracted position of the source (Az/El 3-component vector).
 */
void refv(const vector<double> vu, double refa, double refb, vector<double> vr) {
    // initial estimate = unrefracted vector
    const double x = vu[0];
    const double y = vu[1];
    const double z1 = vu[2];

    // keep correction approximately constant below about 3 deg elevation
    const double z = std::max(z1, 0.05);

    // one Newton-Raphson iteration
    const double z_squared = z * z;
    const double r_squared = x * x + y * y;
    const double r = std::sqrt(r_squared);
    const double wb = refb * r_squared / z_squared;
    const double wt = (refa + wb) / (1.0 + (refa + 3.0 * wb) * (z_squared + r_squared) / z_squared);
    const double d = wt * r / z;
    const double cd = 1.0 - d * d / 2.0;
    const double f = cd * (1.0 - wt);

    // return post-refraction x, y, z
    vr[0] = x * f;
    vr[1] = y * f;
    vr[2] = cd * (z + d * r) + (z1 - z);
}

}
