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
 * Decomposes an [x,y] linear fit into its constituent parameters: zero points, scales, non-perpendicularity and
 * orientation (double precision).
 *
 * The model transforms coordinates [x1,y1] into coordinates [x2,y2] as follows:
 *
 *   x2 = A + B * x1 + C * y1
 *   y2 = D + E * x1 + F * y1
 *
 * The transformation can be decomposed into four steps:
 *
 * 1) Zero points:
 *
 *   x' = xz + x1
 *   y' = yz + y1
 *
 * 2) Scales:
 *
 *   x'' = xs * x'
 *   y'' = ys * y'
 *
 * 3) Non-perpendicularity:
 *
 *   x''' = cos(`perp` / 2) * x'' + sin(`PERP` / 2) * y''
 *   y''' = sin(`perp` / 2) * x'' + cos(`PERP` / 2) * y''
 *
 * 4) Orientation:
 *
 *   x2 = cos(`orient`) * x''' + sin(`orient`) * y'''
 *   y2 = -sin(`orient`) * y''' + cos(`orient`) * y'''
 *
 * See also sla::fitxy(), sla::pxy(), sla::invf(), and sla::xy2xy() functions.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param model Transformation coefficients.
 * @param xz X zero point.
 * @param yz Y zero point.
 * @param xs X scale.
 * @param ys Y scale.
 * @param perp Non-perpendicularity (radians).
 * @param orient Orientation (radians).
 */
void dcmpf(const FitCoeffs& model, double& xz, double& yz, double& xs, double& ys,  double& perp, double& orient) {
    // copy coefficients of the model
    const double a = model.get_a();
    double b = model.get_b();
    const double c = model.get_c();
    const double d = model.get_d();
    double e = model.get_e();
    const double f = model.get_f();

    // scales
    const double rb2e2 = std::sqrt(b * b + e * e);
    const double rc2f2 = std::sqrt(c * c + f * f);
    double xsc;
    if (b * f - c * e >= 0.0) {
        xsc = rb2e2;
    } else {
        b = -b;
        e = -e;
        xsc = -rb2e2;
    }
    const double ysc = rc2f2;

    // non-perpendicularity
    const double p1 = (c != 0.0 || f != 0.0)? std::atan2(c, f): 0.0;
    const double p2 = (e != 0.0 || b != 0.0)? std::atan2(e, b): 0.0;
    const double p = drange(p1 + p2);

    // orientation
    const double ws = c * rb2e2 - e * rc2f2;
    const double wc = b * rc2f2 + f * rb2e2;
    const double orientation = (ws != 0.0 || wc != 0.0)? std::atan2(ws, wc): 0.0;

    // zero points
    const double hp = p / 2.0;
    const double shp = std::sin(hp);
    const double chp = std::cos(hp);
    const double sor = std::sin(orientation);
    const double cor = std::cos(orientation);
    const double det = xsc * ysc * (chp + shp) * (chp - shp);
    double X0, Y0;
    if (std::abs(det) > 0.0) {
        X0 = ysc * (a * (chp * cor - shp * sor) - d * (chp * sor + shp * cor)) / det;
        Y0 = xsc * (a * (chp * sor - shp * cor) + d * (chp * cor + shp * sor)) / det;
    } else {
        X0 = 0.0;
        Y0 = 0.0;
    }

    // return results
    xz = X0;
    yz = Y0;
    xs = xsc;
    ys = ysc;
    perp = p;
    orient = orientation;
}

}
