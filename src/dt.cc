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
 * Estimates and returns the offset between dynamical time and Universal Time for a given historical epoch.
 *
 * Depending on the epoch, one of three parabolic approximations is used:
 *
 *   before 979    Stephenson & Morrison's 390 BC to AD 948 model
 *   979 to 1708   Stephenson & Morrison's 948 to 1600 model
 *   after 1708    McCarthy & Babcock's post-1650 model
 *
 * The breakpoints are chosen to ensure continuity: they occur at places where the adjacent models give the same
 * answer as each other. The models used are based on a lunar tidal acceleration value of -26.00 arcseconds per
 * century.
 *
 * The accuracy is modest, with errors of up to 20 sec during the interval since 1650, rising to perhaps 30 min by
 * 1000 BC. Comparatively accurate values from AD 1600 are tabulated in the Astronomical Almanac (see section K8 of
 * the 1995 AA). The use of double-precision for both argument and result is purely for compatibility with other
 * SLALIB time functions.
 *
 * Reference: Explanatory Supplement to the Astronomical Almanac, ed P.K.Seidelmann, University Science Books (1992),
 *            section 2.553, p83. This contains references to the Stephenson & Morrison and McCarthy & Babcock papers.
 *
 * Original FORTRAN code by Rutherford Appleton Laboratory / P.T. Wallace.
 *
 * @param epoch (Julian) epoch (e.g. 1850.0).
 * @return Rough estimate of ET-UT (after 1984, TT-UT) at the given epoch (seconds).
 */
double dt(double epoch) {
    // centuries since 1800
    const double centuries = (epoch - 1800.0) / 100.0;

    // select model
    if (epoch >= 1708.185161980887) {
        // post-1708: use McCarthy & Babcock
        const double w = centuries - 0.19;
        return 5.156 + 13.3066 * w * w;
    } else if (epoch >= 979.0258204760233) {
        // 979-1708: use Stephenson & Morrison's 948-1600 model
        return 25.5 * centuries * centuries;
    } else {
        // pre-979: use Stephenson & Morrison's 390 BC to AD 948 model
        return 1360.0 + (320.0 + 44.3 * centuries) * centuries;
    }
}

}
