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
 * Auxiliary function used by sla::refro(); calculates refractive index and derivative with respect to height
 * for the stratosphere.
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

/*
 * Auxiliary function used by sla::refro(); calculates refractive index and derivative with respect to height for
 * the troposphere.
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param r0 Height of observer from centre of the Earth (meters).
 * @param t0 Temperature at the observer (degrees K).
 * @param alpha Alpha (see HMNAO paper).
 * @param gamm2 Gamma minus 2 (see HMNAO paper).
 * @param delm2 Delta minus 2 (see HMNAO paper).
 * @param c1 Useful term (see source code of the sla::refro() function).
 * @param c2 Useful term (see source code of the sla::refro() function).
 * @param c3 Useful term (see source code of the sla::refro() function).
 * @param c4 Useful term (see source code of the sla::refro() function).
 * @param c5 Useful term (see source code of the sla::refro() function).
 * @param c6 Useful term (see source code of the sla::refro() function).
 * @param r Current distance from the centre of the Earth (meters).
 * @param t Return value: temperature at `r` (degrees K).
 * @param dn Return value: refractive index at `r`.
 * @param rdndr Return value: `r` * rate the refractive index is changing at `r`.
 */
static void atmt(double r0, double t0, double alpha, double gamm2, double delm2,
    double c1, double c2, double c3, double c4, double c5, double c6, double r,
    double& t, double& dn, double rdndr) {
    t = std::max(std::min(t0 - alpha * (r - r0), 320.0), 100.0);
    const double tt0 = t / t0;
    const double tt0gm2 = std::pow(tt0, gamm2);
    const double tt0dm2 = std::pow(tt0, delm2);
    dn = 1.0 + (c1 * tt0gm2 - (c2 - c5 / t) * tt0dm2) * tt0;
    rdndr = r * (-c3 * tt0gm2 + (c4 - c6 / tt0) * tt0dm2);
}

}
