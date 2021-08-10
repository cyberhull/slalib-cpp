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
 * Adjusts an unrefracted zenith distance to include the effect of atmospheric refraction, using the simple
 * A tan Z + B tan**3 Z model (plus special handling for large ZDs).
 *
 * This function applies the adjustment for refraction in the opposite sense to the usual one - it takes an
 * unrefracted (in vacuo) position and produces an observed (refracted) position, whereas the A tan Z + B tan**3 Z
 * model strictly applies to the case where an observed position is to have the refraction removed. The unrefracted
 * to refracted case is harder, and requires an inverted form of the text-book refraction models; the formula used
 * here is based on the Newton-Raphson method. For the utmost numerical consistency with the refracted to unrefracted
 * model, two iterations are carried out, achieving agreement at the 1D-11 arcseconds level for a ZD of 80 degrees.
 * The inherent accuracy of the model is, of course, far worse than this -- see the documentation for sla::refco()
 * function for more information.
 *
 * At ZD 83 degrees, the rapidly-worsening A tan Z + B tan^3 Z model is abandoned and an empirical formula takes over.
 * For optical/IR wavelengths, over a wide range of observer heights and corresponding temperatures and pressures, the
 * following levels of accuracy (arcseconds, worst case) are achieved, relative to numerical integration through a
 * model atmosphere:
 *
 *   ZR    error
 *   80      0.7
 *   81      1.3
 *   82      2.4
 *   83      4.7
 *   84      6.2
 *   85      6.4
 *   86      8
 *   87     10
 *   88     15
 *   89     30
 *   90     60
 *   91    150      } relevant only to
 *   92    400      } high-elevation sites
 *
 * For radio wavelengths, the errors are typically 50% larger than the optical figures and by ZD 85 degrees are twice
 * as bad, worsening rapidly below that. To maintain 1 arcsecond accuracy down to ZD==85 at the Green Bank site,
 * Condon (2004) has suggested amplifying the amount of refraction predicted by sla::refz() below 10.8 degrees
 * elevation by the factor (1 + 0.00195 * (10.8 - E_t)), where E_t is the unrefracted elevation in degrees.
 *
 * The high-ZD model is scaled to match the normal model at the transition point; there is no glitch.
 *
 * Beyond 93 degrees zenith distance, the refraction is held at its 93 degrees value.
 *
 * See also the sla::refv() function, which performs the adjustment in Cartesian Az/El coordinates, and with the
 * emphasis on speed rather than numerical accuracy.
 *
 * Original FORTRAN code by Rutherford Appleton Laboratory / P.T. Wallace.
 *
 * @param zu Unrefracted zenith distance of the source (radians).
 * @param refa Tan Z coefficient (radians).
 * @param refb Tan**3 Z coefficient (radians).
 * @return Refracted zenith distance (radians).
 */
double refz(double zu, double refa, double refb) {
    // radians to degrees
    constexpr double RADIANS2DEGREES = 57.29577951308232;
    // largest usable ZD (degrees)
    constexpr double DEG93_IN_RADIANS = 93.0;
    // coefficients for high ZD model (used beyond ZD 83 degrees)
    constexpr double C1 = +0.55445;
    constexpr double C2 = -0.01133;
    constexpr double C3 = +0.00202;
    constexpr double C4 = +0.28385;
    constexpr double C5 = +0.02390;
    // ZD at which one model hands over to the other (radians)
    constexpr double ZD_THRESHOLD83 = 83.0 / RADIANS2DEGREES;
    // high-ZD-model prediction (degrees) for that point
    constexpr double REF83 = (C1 + C2 * 7.0 + C3 * 49.0) / (1.0 + C4 * 7.0 + C5 * 49.0);

    // perform calculations for zu or 83 degrees, whichever is smaller
    const double zu1 = std::min(zu, ZD_THRESHOLD83);

    // functions of ZD
    double zl = zu1;
    double sine = std::sin(zl);
    double cosine = std::cos(zl);
    double tangent = sine / cosine;
    double tangent_sqr = tangent * tangent;
    double tangent_cube = tangent * tangent_sqr;

    // refracted ZD (mathematically to better than 1 milliarcsecond at 70 degrees)
    zl = zl - (refa * tangent + refb * tangent_cube) / (1.0 + (refa + 3.0 * refb * tangent_sqr) / (cosine * cosine));

    // further iteration
    sine = std::sin(zl);
    cosine = std::cos(zl);
    tangent = sine / cosine;
    tangent_sqr = tangent * tangent;
    tangent_cube = tangent * tangent_sqr;
    double ref = zu1 - zl + (zl - zu1 + refa * tangent + refb * tangent_cube) /
        (1.0 + (refa + 3.0 * refb * tangent_sqr) / (cosine * cosine));

    // special handling for large zu
    if (zu > zu1) {
        const double E = 90.0 - std::min(DEG93_IN_RADIANS, zu * RADIANS2DEGREES);
        const double E2 = E * E;
        ref = (ref / REF83) * (C1 + C2 * E + C3 * E2) / (1.0 + C4 * E + C5 * E2);
    }
    // return refracted ZD
    return zu - ref;
}

}
