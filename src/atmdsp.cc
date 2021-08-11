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
 * Applies atmospheric-dispersion adjustments to refraction coefficients.
 *
 * To use this routine, first call functio0n sla::refco() specifying `wl1` as the wavelength. This yields refraction
 * coefficients `a1`, `b1`, correct for that wavelength. Subsequently, calls to sla::atmdsp() specifying different
 * wavelengths will produce new, slightly adjusted refraction coefficients which apply to the specified wavelength.
 *
 * Most of the atmospheric dispersion happens between 0.7 micrometre and the UV atmospheric cutoff, and the effect
 * increases strongly towards the UV end. For this reason, a blue reference wavelength is recommended, for example
 * 0.4 micrometres.
 *
 * The accuracy, for this set of conditions:
 *
 *   Height above sea level:  2000 m
 *                 Latitude:  29 deg
 *                 Pressure:  793 mb
 *              Temperature:  17 degC
 *                 Humidity:  50%
 *               Lapse rate:  0.0065 degC/m
 *     Reference wavelength:  0.4 micrometre
 *           Star elevation:  15 deg
 *
 * is about 2.5 mas RMS between 0.3 and 1.0 micrometres, and stays within 4 mas for the whole range longward of 0.3
 * micrometres (compared with a total dispersion from 0.3 to 20.0 micrometres of about 11 arcsec). These errors are
 * typical for ordinary conditions and the given elevation; in extreme conditions, values a few times this size may
 * occur, while at higher elevations the errors become much smaller.
 *
 * If either wavelength exceeds 100 micrometres, the radio case is assumed and the returned refraction coefficients
 * are the same as the given ones. Note that radio refraction coefficients cannot be turned into optical values using
 * this routine, nor vice versa.
 *
 * The algorithm consists of calculation of the refractivity of the air at the observer for the two wavelengths, using
 * the methods of the sla::refro() function, and then scaling of the two refraction coefficients according to classical
 * refraction theory. This amounts to scaling the A coefficient in proportion to (n-1) and the B coefficient almost in
 * the same ratio (see R.M.Green, "Spherical Astronomy", Cambridge University Press, 1985).
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param tdk Ambient temperature (degrees K).
 * @param pmb Ambient pressure (millibars).
 * @param rh Ambient relative humidity (0.0..1.0).
 * @param wl1 Reference wavelength (micrometres); 0.4D0 recommended.
 * @param a1 Refraction coefficient A for wavelength `wl1` (radians).
 * @param b1 Refraction coefficient B for wavelength `wl1` (radians).
 * @param wl2 Wavelength for which adjusted A, B required.
 * @param a2 Return value: refraction coefficient A for wavelength `wl2` (radians).
 * @param b2 Return value: refraction coefficient B for wavelength `wl2` (radians).
 */
void atmdsp(double tdk, double pmb, double rh, double wl1, double a1, double b1, double wl2, double& a2, double& b2) {
    // check for radio wavelengths
    if (wl1 > 100.0 || wl2 > 100.0) {
        // radio: no dispersion
        a2 = a1;
        b2 = b1;
    } else {
        // optical: keep arguments within safe bounds
        const double tdkok = std::min(std::max(tdk, 100.0), 500.0);
        const double pmbok = std::min(std::max(pmb, 0.0), 10000.0);
        const double rhok = std::min(std::max(rh, 0.0), 1.0);

        // atmosphere parameters at the observer
        const double psat = std::pow(10.0, (-8.7115 + 0.03477 * tdkok));
        const double pwo = rhok * psat;
        const double w1 = 11.2684e-6 * pwo;

        // refractivity at the observer for first wavelength
        double wlok = std::max(wl1, 0.1);
        double wlsq = wlok * wlok;
        double w2 = 77.5317e-6 + (0.43909e-6 + 0.00367e-6 / wlsq) / wlsq;
        const double dn1 = (w2 * pmbok - w1) / tdkok;

        // refractivity at the observer for second wavelength
        wlok = std::max(wl2, 0.1);
        wlsq = wlok * wlok;
        w2 = 77.5317e-6 + (0.43909e-6 + 0.00367e-6 / wlsq) / wlsq;
        const double dn2 = (w2 * pmbok - w1) / tdkok;

        // scale the refraction coefficients (see Green 4.31, p93)
        if (dn1 != 0.0) {
            const double f = dn2 / dn1;
            a2 = a1 * f;
            b2 = b1 * f;
            if (dn1 != a1) {
                b2 *= (1.0 + dn1 * (dn1 - dn2) / (2.0 * (dn1 - a1)));
            }
        } else {
            a2 = a1;
            b2 = b1;
        }
    }
}

}
