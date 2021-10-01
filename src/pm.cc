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
 * Applies corrections for proper motion to a star RA,Dec (double precision).
 *
 * If the available proper motions are pre-FK5 they will be per tropical year rather than per Julian year, and so the
 * epochs must both be Besselian rather than Julian; in such cases, a scaling factor of 365.2422D0/365.25D0 should be
 * applied to the radial velocity before use.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param dir_ep0 RA,Dec at epoch EP0 (radians).
 * @param motion Proper motions: RA,Dec changes per year of epoch; the proper motions in RA are dRA/dt rather than
 *   cos(Dec)*dRA/dt, and are in the same coordinate system as `dir_ep0`.
 * @param parallax Parallax (arcseconds).
 * @param r_velocity Radial velocity (km/sec, positive if receding).
 * @param ep0 Start epoch in years (e.g. Julian epoch).
 * @param ep1 End epoch in years (same system as `ep0`).
 * @param dir_ep1 Return value: RA,Dec at epoch `ep1` (radians).
 */
void pm(const Spherical<double>& dir_ep0, const Spherical<double>& motion, double parallax, double r_velocity,
    double ep0, double ep1, Spherical<double>& dir_ep1) {
    // km/s to AU/year multiplied by arcseconds to radians
    constexpr double VFR = (365.25 * 86400.0 / 149597870.0) * 4.8481368111e-6;

    // spherical to Cartesian
    Vector<double> pos;
    dcs2c(dir_ep0, pos);

    // space motion (radians per year)
    const double w = VFR * r_velocity * parallax;
    const Vector<double> em = {
        -motion.get_ra() * pos[1] - motion.get_dec() * std::cos(dir_ep0.get_ra()) * std::sin(dir_ep0.get_dec()) + w * pos[0],
        motion.get_ra() * pos[0] - motion.get_dec() * std::sin(dir_ep0.get_ra()) * std::sin(dir_ep0.get_dec()) + w * pos[1],
        motion.get_dec() * std::cos(dir_ep0.get_dec()) + w * pos[2]
    };

    // apply the motion
    const double time = ep1 - ep0;
    for (int i = 0; i < 3; i++) {
        pos[i] += time * em[i];
    }

    // Cartesian to spherical
    dcc2s(pos, dir_ep1);
    dir_ep1.set_ra(dranrm(dir_ep1.get_ra()));
}

}
