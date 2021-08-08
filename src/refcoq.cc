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
 * Determines the constants A and B in the atmospheric refraction model dZ = A tan Z + B tan**3 Z. This is a fast
 * alternative to the sla::refco() function.
 *
 * Z is the "observed" zenith distance (i.e. affected by refraction) and dZ is what to add to Z to give the
 * "topocentric" (i.e. in vacuo) zenith distance.
 *
 * The model is an approximation, for moderate zenith distances, to the predictions of the sla::refro() function. The
 * approximation is maintained across a range of conditions, and applies to both optical/IR and radio.
 *
 * The algorithm is a fast alternative to the sla::refco function. The latter calls the sla::refro() function itself:
 * this involves integrations through a model atmosphere, and is costly in processor time. However, the model which is
 * produced is precisely correct for two zenith distance (45 degrees and about 76 degrees) and at other zenith
 * distances is limited in accuracy only by the A tan Z + B tan**3 Z formulation itself. The present routine is not as
 * accurate, though it satisfies most practical requirements.
 *
 * The model omits the effects of (i) height above sea level (apart from the reduced pressure itself), (ii) latitude
 * (i.e. the flattening of the Earth) and (iii) variations in tropospheric lapse rate.
 *
 *  The model was tested using the following range of conditions:
 *
 *    Lapse rates 0.0055, 0.0065, 0.0075 K/metre
 *    Latitudes 0, 25, 50, 75 degrees
 *    Heights 0, 2500, 5000 metres ASL
 *    Pressures mean for height -10% to +5% in steps of 5%
 *    Temperatures -10 deg to +20 deg with respect to 280 deg at SL
 *    Relative humidity 0, 0.5, 1
 *    Wavelengths 0.4, 0.6, ... 2 micron, + radio
 *    Zenith distances 15, 45, 75 degrees
 *
 * The accuracy with respect to direct use of the sla_REFRO routine was as follows:
 *
 *               Worst    RMS
 *   Optical/IR  62 mas   8 mas
 *   Radio       319 mas  49 mas
 *
 * For this particular set of conditions:
 *
 *   Lapse rate:   0.0065 K/metre
 *   Latitude:     50 degrees
 *   Altitude:     sea level
 *   Pressure:     1005 mb
 *   Temperature:  280.15 K
 *   Humidity:     80%
 *   Wavelength:   5740 Angstroms
 *
 * the results were as follows:
 *
 *   ZD   sla_REFRO  sla_REFCOQ  Saastamoinen
 *   deg  arcsec     arcsec      arcsec
 *   10    10.27      10.27       10.27
 *   20    21.19      21.20       21.19
 *   30    33.61      33.61       33.60
 *   40    48.82      48.83       48.81
 *   45    58.16      58.18       58.16
 *   50    69.28      69.30       69.27
 *   55    82.97      82.99       82.95
 *   60   100.51     100.54      100.50
 *   65   124.23     124.26      124.20
 *   70   158.63     158.68      158.61
 *   72   177.32     177.37      177.31
 *   74   200.35     200.38      200.32
 *   76   229.45     229.43      229.42
 *   78   267.44     267.29      267.41
 *   80   319.13     318.55      319.10
 *
 * The values for Saastamoinen's formula (which includes terms up to tan^5) are taken from Hohenkerk & Sinclair (1985).
 *
 * The results from the much slower but more accurate sla::refco() function have not been included in the tabulation
 * as they are identical to those in the sla::refro() column to the 0.01 arcsec resolution used.
 *
 * The algorithm draws on several sources, as follows:
 *
 * a) The formula for the saturation vapour pressure of water as a function of temperature and temperature is taken
 *   from expressions A4.5-A4.7 of Gill (1982).
 *
 * b) The formula for the water vapour pressure, given the saturation pressure and the relative humidity, is from
 *   Crane (1976), expression 2.5.5.
 *
 * c) The refractivity of air is a function of temperature, total pressure, water-vapour pressure and, in the case
 *   of optical/IR but not radio, wavelength. The formulae for the two cases are developed from Hohenkerk & Sinclair
 *   (1985) and Rueger (2002).
 *
 * The above three items are as used in the sla::refro() function.
 *
 * d) The formula for beta, the ratio of the scale height of the atmosphere to the geocentric distance of the
 *   observer, is an adaption of expression 9 from Stone (1996). The adaptations, arrived at empirically, consist
 *   of (i) a small adjustment to the coefficient and (ii) a humidity term for the radio case only.
 *
 * e) The formulae for the refraction constants as a function of n-1 and beta are from Green (1987), expression 4.31.
 *
 * References:
 *
 *   Crane, R.K., Meeks, M.L. (ed), "Refraction Effects in the Neutral Atmosphere", Methods of Experimental Physics:
 *   Astrophysics 12B, Academic Press, 1976.
 *
 *   Gill, Adrian E., "Atmosphere-Ocean Dynamics", Academic Press, 1982.
 *
 *   Green, R.M., "Spherical Astronomy", Cambridge University Press, 1987.
 *
 *   Hohenkerk, C.Y., & Sinclair, A.T., NAO Technical Note No. 63, 1985.
 *
 *   Rueger, J.M., "Refractive Index Formulae for Electronic Distance Measurement with Radio and Millimetre Waves",
 *   in Unisurv Report S-68, School of Surveying and Spatial Information Systems, University of New South Wales,
 *   Sydney, Australia, 2002.
 *
 *   Stone, Ronald C., P.A.S.P. 108 1051-1058, 1996.
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param tdk Ambient temperature at the observer (degrees K); clamped to the [100K..500K] range.
 * @param pmb Pressure at the observer (millibars); clamped to the [0..10000] range; zero pressure yields zero results.
 * @param rh Relative humidity at the observer (range 0..1); clamped to the [0..1] range.
 * @param wl Effective wavelength of the source (micrometres); clamped internally to the [0.1..1E6] range.
 * @param refa Return value: tan Z coefficient (radians).
 * @param refb Return value: tan**3 Z coefficient (radians).
 */
void refcoq(double tdk, double pmb, double rh, double wl, double& refa, double& refb) {
    // decide whether optical/IR or radio case: switch at 100 microns
    const bool optic = wl <= 100.0;

    // restrict parameters to safe values
    const double t = std::min(std::max(tdk, 100.0), 500.0);
    const double p = std::min(std::max(pmb, 0.0), 10000.0);
    const double r = std::min(std::max(rh, 0.0), 1.0);
    const double w = std::min(std::max(wl, 0.1), 1.0e6);

    // water vapour pressure at the observer
    double pw;
    if (p > 0.0) {
        const double tdc = t - 273.15;
        const double ps = std::pow(10.0, ((0.7859 + 0.03477 * tdc) / (1.0 + 0.00412 * tdc))) *
            (1.0 + p * (4.5 - 6.0 + 6.0 - 10.0 * tdc * tdc));
        pw = r * ps / (1.0 - (1.0 - r) * ps / p);
    } else {
        pw = 0.0;
    }
    // refractive index minus 1 at the observer
    double gamma;
    if (optic) {
        const double wlsq = w * w;
        gamma = ((77.53484 - 6.0 + (4.39108 - 7.0 + 3.666 - 9.0 / wlsq) / wlsq) * p - 11.2684 - 6.0 * pw) / t;
    } else {
        gamma = (77.6890 - 6.0 * p - (6.3938 - 6.0 - 0.375463 / t) * pw) / t;
    }
    // formula for beta adapted from Stone, with empirical adjustments
    double beta = 4.4474 - 6.0 * t;
    if (!optic) {
        beta -= 0.0074 * pw * beta;
    }
    // refraction constants from Green
    refa = gamma * (1.0 - beta);
    refb = -gamma * (beta - gamma / 2.0);
}

}
