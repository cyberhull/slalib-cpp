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
static void atms(double rt, double tt, double dnt, double gamal, double r, double& dn, double& rdndr) {
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
 * @param c5 Useful term; zero in the optical case (see source code of the sla::refro() function).
 * @param c6 Useful term; zero in the optical case (see source code of the sla::refro() function).
 * @param r Current distance from the centre of the Earth (meters).
 * @param t Return value: temperature at `r` (degrees K).
 * @param dn Return value: refractive index at `r`.
 * @param rdndr Return value: `r` * rate the refractive index is changing at `r`.
 */
static void atmt(double r0, double t0, double alpha, double gamm2, double delm2,
    double c1, double c2, double c3, double c4, double c5, double c6, double r,
    double& t, double& dn, double& rdndr) {
    t = std::max(std::min(t0 - alpha * (r - r0), 320.0), 100.0);
    const double tt0 = t / t0;
    const double tt0gm2 = std::pow(tt0, gamm2);
    const double tt0dm2 = std::pow(tt0, delm2);
    dn = 1.0 + (c1 * tt0gm2 - (c2 - c5 / t) * tt0dm2) * tt0;
    rdndr = r * (-c3 * tt0gm2 + (c4 - c6 / tt0) * tt0dm2);
}

/**
 * Calculates atmospheric refraction for radio and optical/IR wavelengths.
 *
 * This function computes the refraction for zenith distances up to and a little beyond 90 degrees using the method
 * of Hohenkerk and Sinclair (NAO Technical Notes 59 and 63, subsequently adopted in the Explanatory Supplement, 1992
 * edition - see section 3.281).
 *
 * As in the original Hohenkerk and Sinclair algorithm, fixed values of the water vapour polytrope exponent, the height
 * of the tropopause, and the height at which refraction is negligible are used.
 *
 * The radio refraction has been tested against work done by Iain Coulson, JACH, (private communication 1995) for the
 * James Clerk Maxwell Telescope, Mauna Kea. For typical conditions, agreement at the 0.1 arcsec level is achieved for
 * moderate ZD, worsening to perhaps 0.5-1.0 arcsec at ZD 80 deg. At hot and humid sea-level sites the accuracy will
 * not be as good.
 *
 * The algorithm is designed for observers in the troposphere. The supplied temperature, pressure and lapse rate are
 * assumed to be for a point in the troposphere and are used to define a model atmosphere with the tropopause at 11km
 * altitude and a constant temperature above that. However, in practice, the refraction values returned for
 * stratospheric observers, at altitudes up to 25km, are quite usable.
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * The FORTRAN code was a development of the optical/IR refraction subroutine AREF of C.Hohenkerk (HMNAO, September
 * 1984), with extensions to support the radio case. Apart from merely cosmetic changes, the following modifications
 * to the original HMNAO optical/IR refraction code had been made:
 *
 * - The angle arguments have been changed to radians.
 *
 * - Any value of ZOBS is allowed.
 *
 * - Other argument values have been limited to safe values.
 *
 * - Murray's values for the gas constants have been used (Vectorial Astrometry, Adam Hilger, 1983).
 *
 * - The numerical integration phase has been rearranged for extra clarity.
 *
 * - A better model for Ps(T) has been adopted (taken from Gill, Atmosphere-Ocean Dynamics, Academic Press, 1982).
 *
 * - More accurate expressions for Pwo have been adopted (again from Gill 1982).
 *
 * - The formula for the water vapour pressure, given the saturation pressure and the relative humidity, is from
 *   Crane (1976), expression 2.5.5.
 *
 * - Provision for radio wavelengths has been added using expressions devised by A.T.Sinclair, RGO (private
 *   communication 1989). The refractivity model currently used is from J.M.Rueger, "Refractive Index Formulae for
 *   Electronic Distance Measurement with Radio and Millimetre Waves", in Unisurv Report S-68 (2002), School of
 *   Surveying and Spatial Information Systems, University of New South Wales, Sydney, Australia.
 *
 * - The optical refractivity for dry air is from Resolution 3 of the International Association of Geodesy adopted at
 *   the XXIIth General Assembly in Birmingham, UK, 1999.
 *
 * - Various small changes have been made to gain speed.
 *
 * @param zobs Observed zenith distance of the source (radians); before use, the value of `zobs` is expressed in the
 *   range +/- pi; if this ranged `zobs` is -ve, the return value is computed from its absolute value before being
 *   made -ve to match; in addition, if it has an absolute value greater than 93 degrees, a fixed refraction value
 *   equal to the result for `zobs` = 93 degrees is returned, appropriately signed.
 * @param hm Height of the observer above sea level (meters).
 * @param tdk Ambient temperature at the observer (degrees K).
 * @param pmb Pressure at the observer (millibars).
 * @param rh Relative humidity at the observer (range: [0-1]); relative humidity `rh` is formally defined in terms of
 *   "mixing ratio" rather than pressures or densities as is often stated. It is the mass of water per unit mass of
 *   dry air divided by that for saturated air at the same temperature and pressure (see Gill 1982).
 * @param wl Effective wavelength of the source (micrometers); the radio refraction is chosen by specifying WL > 100
 *   micrometers; because the algorithm takes no account of the ionosphere, the accuracy deteriorates at low
 *   frequencies, below about 30 MHz.
 * @param phi Latitude of the observer (radians, astronomical).
 * @param tlr Temperature lapse rate in the troposphere (degrees K/meter); A suggested value for the TLR argument
 *   is 0.0065; the refraction is significantly affected by `tlr`, and if studies of the local atmosphere have been
 *   carried out a better `tlr` value may be available; the sign of the supplied TLR value is ignored.
 * @param eps Precision required to terminate iteration (radians); a suggested value for the EPS argument is 1.0e-8;
 *   the result is usually at least two orders of magnitude more computationally precise than the supplied `eps` value.
 * @return Refraction: in vacuo ZD minus observed ZD (radians).
 */
double refro(double zobs, double hm, double tdk, double pmb, double rh, double wl, double phi, double tlr, double eps) {
    // 93 degrees in radians
    constexpr double DEG93_IN_RADIANS = 1.623156204;
    // universal (molar) gas constant
    constexpr double MOLAR_GAS_CONSTANT = 8314.32;
    // molecular weight of dry air
    constexpr double DRY_AIR_MOL_WEIGHT = 28.9644;
    // molecular weight of water vapour
    constexpr double WATER_VAPOUR_MOL_WEIGHT = 18.0152;
    // mean Earth radius (meters)
    constexpr double EARTH_RADIUS = 6378120.0;
    // exponent of temperature dependence of water vapour pressure
    constexpr double DELTA = 18.36;
    // height of tropopause (meters)
    constexpr double TROPOPAUSE_HEIGHT = 11000.0;
    // upper limit for refractive effects (meters)
    constexpr double RE_HEIGHT_LIMIT = 80000.0;
    // numerical integration: maximum number of strips
    constexpr int MAX_STRIPS = 16384;

    // the refraction integrand
    auto refraction_integrand = [] (double dn, double rdndr) -> double {
        return rdndr / (dn + rdndr);
    };

    // transform zobs into the normal range
    const double zobs1 = drange(zobs);
    const double zobs2 = std::min(std::abs(zobs1), DEG93_IN_RADIANS);

    // keep other arguments within safe bounds
    const double hm_ok = std::min(std::max(hm, -1.0e3), RE_HEIGHT_LIMIT);
    const double tdk_ok = std::min(std::max(tdk, 100.0), 500.0);
    const double pmb_ok = std::min(std::max(pmb, 0.0), 10000.0);
    const double rh_ok = std::min(std::max(rh, 0.0), 1.0);
    const double wl_ok = std::max(wl, 0.1);
    const double alpha = std::min(std::max(std::abs(tlr), 0.001), 0.01);

    // tolerance for iteration
    const double tolerance = std::min(std::max(std::abs(eps), 1.0e-12), 0.1) / 2.0;

    // decide whether optical/IR or radio case - switch at 100 microns
    bool optic = wl_ok <= 100.0;

    // set up model atmosphere parameters defined at the observer
    const double wl_squared = wl_ok * wl_ok;
    const double gb = 9.784 * (1.0 - 0.0026 * std::cos(phi + phi) - 0.00000028 * hm_ok);
    const double a = optic?
        (287.6155 + (1.62887 + 0.01360 / wl_squared) / wl_squared) * 273.15e-6 / 1013.25: 77.6890e-6;
    const double gamal = (gb * DRY_AIR_MOL_WEIGHT) / MOLAR_GAS_CONSTANT;
    const double gamma = gamal / alpha;
    const double gamm2 = gamma - 2.0;
    const double delm2 = DELTA - 2.0;
    const double tdc = tdk_ok - 273.15;
    const double psat = std::pow(10.0, ((0.7859 + 0.03477 * tdc) / (1.0 + 0.00412 * tdc))) *
        (1.0 + pmb_ok * (4.5e-6 + 6.0e-10 * tdc * tdc));
    const double pwo = pmb_ok > 0.0 ? rh_ok * psat / (1.0 - (1.0 - rh_ok) * psat / pmb_ok) : 0.0;
    const double w = pwo * (1.0 - WATER_VAPOUR_MOL_WEIGHT / DRY_AIR_MOL_WEIGHT) * gamma / (DELTA - gamma);
    const double c1 = a * (pmb_ok + w) / tdk_ok;
    const double c2 = (a * w + (optic ? 11.2684e-6 : 6.3938e-6) * pwo) / tdk_ok;
    const double c3 = (gamma - 1.0) * alpha * c1 / tdk_ok;
    const double c4 = (DELTA - 1.0) * alpha * c2 / tdk_ok;
    const double c5 = optic ? 0.0 : 375463.0e-6 * pwo / tdk_ok;
    const double c6 = optic ? 0.0 : c5 * delm2 * alpha / (tdk_ok * tdk_ok);

    // conditions at the observer
    const double r0 = EARTH_RADIUS + hm_ok;
    double temp0, dn0, rdndr0;
    atmt(r0, tdk_ok, alpha, gamm2, delm2, c1, c2, c3, c4, c5, c6, r0, temp0, dn0, rdndr0);
    const double sk0 = dn0 * r0 * std::sin(zobs2);
    const double f0 = refraction_integrand(dn0, rdndr0);

    // conditions in the troposphere at the tropopause
    const double rt = EARTH_RADIUS + std::max(TROPOPAUSE_HEIGHT, hm_ok);
    double tt, dnt, rdndrt;
    atmt(r0, tdk_ok, alpha, gamm2, delm2, c1, c2, c3, c4, c5, c6, rt, tt, dnt, rdndrt);
    double sine = sk0 / (rt * dnt);
    const double zt = std::atan2(sine, std::sqrt(std::max(1.0 - sine * sine, 0.0)));
    const double ft = refraction_integrand(dnt, rdndrt);

    // conditions in the stratosphere at the tropopause
    double dnts, rdndrp;
    atms(rt, tt, dnt, gamal, rt, dnts, rdndrp);
    sine = sk0 / (rt * dnts);
    const double zts = std::atan2(sine, std::sqrt(std::max(1.0 - sine * sine, 0.0)));
    const double fts = refraction_integrand(dnts, rdndrp);

    // conditions at the stratosphere limit
    const double rs = EARTH_RADIUS + RE_HEIGHT_LIMIT;
    double dns, rdndrs;
    atms(rt, tt, dnt, gamal, rs, dns, rdndrs);
    sine = sk0 / (rs * dns);
    const double zs = std::atan2(sine, std::sqrt(std::max(1.0 - sine * sine, 0.0)));
    const double fs = refraction_integrand(dns, rdndrs);

    // variable initialization to avoid compiler warning
    double refp, reft = 0.0;

    // integrate the refraction integral in two parts; first in the troposphere (K=1), then in the stratosphere (K=2)
    for (int k = 0; k < 2; k++) {
        // initialize previous refraction to ensure at least two iterations
        double ref_old = 1.0;
        // start off with 8 strips
        int num_strips = 8;
        // start Z, Z range, and start and end values
        double z0, z_range, fb, ff;
        if (k == 0) {
            z0 = zobs2;
            z_range = zt - z0;
            fb = f0;
            ff = ft;
        } else {
            z0 = zts;
            z_range = zs - z0;
            fb = fts;
            ff = fs;
        }
        // sums of odd and even values
        double f_odd = 0.0;
        double f_even = 0.0;

        // start of iteration loop (terminates at specified precision)
        for (int step = 1;;) { // first time through the loop we have to do every point
            // strip width
            const double h = z_range / (double) num_strips;
            // initialize distance from Earth centre for quadrature pass
            double r = (k == 0) ? r0 : rt;

            // one pass (no need to compute evens after first time)
            for (int i = 1; i <= num_strips - 1; i += step) {
                // sine of observed zenith distance.
                const double sine_zd = std::sin(z0 + h * (double)(i));

                // find r (to the nearest meter, maximum four iterations)
                double dn, rdndr;
                if (sine_zd > 1.0e-20) {
                    const double ww = sk0 / sine_zd;
                    double rg = r;
                    double dr = 1.0e6;
                    for (int j = 0; j < 4 && std::abs(dr) > 1.0; j++) {
                        if (k == 0) {
                            double tg;
                            atmt(r0, tdk_ok, alpha, gamm2, delm2, c1, c2, c3, c4, c5, c6, rg, tg, dn, rdndr);
                        } else {
                            atms(rt, tt, dnt, gamal, rg, dn, rdndr);
                        }
                        dr = (rg * dn - ww) / (dn + rdndr);
                        rg = rg - dr;
                    }
                    r = rg;
                }

                // find the refractive index and integrand at r
                if (k == 0) {
                    double t;
                    atmt(r0, tdk_ok, alpha, gamm2, delm2, c1, c2, c3, c4, c5, c6, r, t, dn, rdndr);
                } else {
                    atms(rt, tt, dnt, gamal, r, dn, rdndr);
                }
                const double f = refraction_integrand(dn, rdndr);

                // accumulate odd and (first time only) even values
                if (step == 1 && i % 2 == 0) {
                    f_even = f_even + f;
                } else {
                    f_odd = f_odd + f;
                }
            }

            // evaluate the integrand using Simpson's Rule
            refp = h * (fb + 4.0 * f_odd + 2.0 * f_even + ff) / 3.0;

            // has the required precision been achieved (or can't be)?
            if (std::abs(refp - ref_old) > tolerance && num_strips < MAX_STRIPS) {
                // NO: prepare for next iteration
                // save current value for convergence test
                ref_old = refp;
                // double the number of strips
                num_strips += num_strips;
                // sum of all current values = sum of next pass's even values
                f_even = f_even + f_odd;
                // prepare for new odd values
                f_odd = 0.0;
                // skip even values next time
                step = 2;
            } else {
                // YES: save troposphere component and terminate the loop
                if (k == 0) {
                    reft = refp;
                }
                break;
            }
        }
    }
    // return result
    double result = reft + refp;
    if (zobs1 < 0.0) {
        result = -result;
    }
    return result;
}

}
