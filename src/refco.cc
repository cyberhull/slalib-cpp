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
 * Determines the constants A and B in the atmospheric refraction model dZ = A tan Z + B tan**3 Z.
 *
 * Z is the "observed" zenith distance (i.e. affected by refraction) and dZ is what to add to Z to give the
 * "topocentric" (i.e. in vacuo) zenith distance.
 *
 * This function is a slower but more accurate alternative to the sla::refcoq() function. The constants it produces
 * give perfect agreement with sla::refro() at zenith distances atan(1.0) (45 degrees) and atan(4.0) (about 76
 * degrees). It achieves 0.5 arcsecond accuracy for ZD < 80 degrees, 0.01 arcsecond accuracy for ZD < 60 deg, and
 * 0.001 arcsecond accuracy for ZD < 45 degrees.
 *
 * Original FORTRAN code by Rutherford Appleton Laboratory / P.T. Wallace.
 *
 * @param hm Height of the observer above sea level (meters).
 * @param tdk Ambient temperature at the observer (degrees K).
 * @param pmb Pressure at the observer (millibars).
 * @param rh Relative humidity at the observer (range: [0..1])
 * @param wl Effective wavelength of the source (micrometers).
 * @param phi Latitude of the observer (radians, astronomical).
 * @param tlr Temperature lapse rate in the troposphere (degrees K/meter).
 * @param eps Precision required to terminate iteration (radians).
 * @param refa Return value: tan Z coefficient (radians).
 * @param refb Return value: tan**3 Z coefficient (radians).
 */
void refco(double hm, double tdk, double pmb, double rh, double wl, double phi, double tlr, double eps,
    double& refa, double& refb) {

    // sample zenith distances: atan(1.0) and atan(4.0)
    constexpr double ATAN_1 = 0.7853981633974483;
    constexpr double ATAN_4 = 1.325817663668033;

    // determine refraction for the two sample zenith distances
    const double r1 = refro(ATAN_1, hm, tdk, pmb, rh, wl, phi, tlr, eps);
    const double r2 = refro(ATAN_4, hm, tdk, pmb, rh, wl, phi, tlr, eps);

    // solve for refraction constants
    refa = (64.0 * r1 - r2) / 60.0;
    refb = (r2 - 4.0 * r1) / 60.0;
}

}
