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

/*
 *  Auxiliary function used by sla::refro(); calculates refractive index and derivative with respect to height
 *  for the stratosphere.
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param rt Height of tropopause from centre of the Earth (meters).
 * @param tt Temperature at the tropopause (degrees K).
 * @param dnt Refractive index at the tropopause.
 * @param gamal Constant of the atmospheric model = G * MD / R.
 * @param r Current distance from the centre of the Earth (meters).
 * @param dn Return value: refractive index at `r`.
 * @param rdndr Return value: `r` * rate the refractive index is changing at `r`.
 */
static void atms(double rt, double tt, double dnt, double gamal, double r, double& dn, double rdndr) {
    const double b = gamal / tt;
    const double w = (dnt - 1.0) * std::exp(-b * (r - rt));
    dn = 1.0 + w;
    rdndr = -r * b * w;
}

}
