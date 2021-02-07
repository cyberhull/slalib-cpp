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
 * Function that calculates air mass at given zenith distance.
 *
 * Uses Hardie's (1962) polynomial fit to Bemporad's data for the relative air mass, X, in units of thickness at the
 * zenith as tabulated by Schoenberg (1929). This is adequate for all normal needs as it is accurate to better than
 * 0.1% up to X = 6.8 and better than 1% up to X = 10. Bemporad's tabulated values are unlikely to be trustworthy to
 * such accuracy because of variations in density, pressure and other conditions in the atmosphere from those assumed
 * in his work.
 *
 * References:
 *   Hardie, R.H., 1962, in "Astronomical Techniques" ed. W.A. Hiltner, University of Chicago Press, p180.
 *   Schoenberg, E., 1929, Hdb. d. Ap., Berlin, Julius Springer, 2, 268.
 *
 *  Original FORTRAN code by P.W.Hill, St Andrews, P.T.Wallace
 *
 * @param zenith_dist Observed zenith distance (radians). Here, "observed" means "as affected by refraction". The sign
 *   of the distance is ignored. At zenith distances greater than about 87 degrees the air mass is held constant to
 *   avoid arithmetic overflows.
 *
 * @return An estimate of the air mass at given zenith distance, in units of that at the zenith.
 */
double airmas(double zenith_dist) {
  const double seczm1 = 1.0 / (std::cos(std::fmin(1.52, std::fabs(zenith_dist)))) - 1.0;
  return 1.0 + seczm1 * (0.9981833 - seczm1 * (0.002875 + 0.0008083 * seczm1));
}

}
