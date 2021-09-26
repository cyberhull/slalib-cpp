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
 * Inverts a linear model of the type produced by the sla::fitxy() function.
 *
 * The coefficients (accessible as methods `get_a()`, `get_b()`, etc.) of the models relate two sets of [x,y]
 * coordinates ([x1,y1] and [x2,y1]) as follows:
 *
 *   x2 = A + B * x1 + C * y1
 *   y2 = D + E * x1 + F * y1
 *
 * The `sla::invf()` function generates a new set of coefficients (A', B', etc.) so that:
 *
 *   x1 = A' + B' * x2 + C' * y2
 *   y1 = D' + E' * x2 + F' * y2
 *
 * Two successive calls to sla::invf() will thus deliver a set of coefficients equal to the starting values.
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param model Model coefficients.
 * @param inverse Inverse model.
 * @return `true` if successful, `false` if there is no inverse transform.
 */
bool invf(const FitCoeffs& model, FitCoeffs& inverse) {
    const double b = model.get_b();
    const double c = model.get_c();
    const double e = model.get_e();
    const double f = model.get_f();
    const double det = b * f - c * e;
    if (det != 0.0) {
        const double a = model.get_a();
        const double d = model.get_d();
        inverse[0] = (c * d - a * f) / det;
        inverse[1] = f / det;
        inverse[2] = -c / det;
        inverse[3] = (a * e - b * d) / det;
        inverse[4] = -e / det;
        inverse[5] = b / det;
        return true;
    } else {
        return false;
    }
}

}
