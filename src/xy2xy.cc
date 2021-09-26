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
 * Transforms one [x,y] into another using a linear model of the type produced by the sla::fitxy() routine (double
 * precision).
 *
 * The model relates two sets of [x,y] coordinates using A, B, C, D, and F coefficients stored in `model`
 * transformation coefficients as follows:
 *
 *   `x2` = A + B * `x1` + C * `y1`
 *   `y2` = D + E * `x1` + F * `y1`
 *
 * See also sla::fitxy(), sla::pxy(), sla::invf(), and sla::dcmpf() functions.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param x1 X coordinate.
 * @param y1 Y coordinate.
 * @param model Transformation coefficients.
 * @param x2 Return value: transformed X coordinate.
 * @param y2 Return value: transformed Y coordinate.
 */
void xy2xy(double x1, double y1, const FitCoeffs& model, double& x2, double& y2) {
    x2 = model.get_a() + model.get_b() * x1 + model.get_c() * y1;
    y2 = model.get_d() + model.get_e() * x1 + model.get_f() * y1;
}

}
