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
 * Forms the matrix of precession between two epochs, using the model of Simon et al (1994), which is suitable for
 * long periods of time (double precision).
 *
 *  The absolute accuracy of the model is limited by the uncertainty in the general precession, about 0.3 arcseconds
 *  per 1000 years. The remainder of the formulation provides a precision of 1 milliarcsecond over the interval from
 *  1000AD to 3000AD, 0.1 arcseconds from 1000BC to 5000AD, and 1 arcsecond from 4000BC to 8000AD.
 *
 * Reference:
 *   Simon, J.L. et al., 1994. Astron.Astrophys., 282, 663-683.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param ep0 Beginning epoch: TDB (Barycentric Dynamical Time) Julian epoch.
 * @param ep1 Ending epoch: TDB (Barycentric Dynamical Time) Julian epoch.
 * @param mat Return value: precession matrix; the matrix is in the sense v(ep1) == mat * v(ep0).
 */
void precl(double ep0, double ep1, Matrix<double> mat) {
    // arc seconds to radians
    constexpr double ARCSECS_2_RADIANS = 0.484813681109535994e-5;

    // interval between basic epoch J2000.0 and beginning epoch (1000JY)
    const double t0 = (ep0 - 2000.0) / 1000.0;

    // interval over which precession required (1000JY)
    const double t = (ep1 - ep0) / 1000.0;

    // euler angles
    const double tas2r = t * ARCSECS_2_RADIANS;
    const double w = 23060.9097 +
                      (139.7459 +
                       (-0.0038 +
                       (-0.5918 +
                       (-0.0037 +
                         0.0007 * t0) * t0) * t0) * t0) * t0;

    const double zeta = (w + (30.2226 +
                             (-0.2523 +
                             (-0.3840 +
                             (-0.0014 +
                               0.0007 * t0) * t0) * t0) * t0 +
                             (18.0183 +
                             (-0.1326 +
                              (0.0006 +
                               0.0005 * t0) * t0) * t0 +
                             (-0.0583 +
                             (-0.0001 +
                               0.0007 * t0) * t0 +
                             (-0.0285 +
                              -0.0002 * t) * t) * t) * t) * t) * tas2r;

    const double z = (w + (109.5270 +
                            (0.2446 +
                           (-1.3913 +
                           (-0.0134 +
                             0.0026 * t0) * t0) * t0) * t0 +
                           (18.2667 +
                           (-1.1400 +
                           (-0.0173 +
                             0.0044 * t0) * t0) * t0 +
                           (-0.2821 +
                           (-0.0093 +
                             0.0032 * t0) * t0 +
                           (-0.0301 +
                             0.0006 * t0
                            -0.0001 * t) * t) * t) * t) * t) * tas2r;

    const double theta = (20042.0207 +
                           (-85.3131 +
                            (-0.2111 +
                             (0.3642 +
                             (0.0008 +
                             -0.0005 * t0) * t0) * t0) * t0) * t0 +
                           (-42.6566 +
                            (-0.2111 +
                             (0.5463 +
                             (0.0017 +
                             -0.0012 * t0) * t0) * t0) * t0 +
                           (-41.8238 +
                             (0.0359 +
                             (0.0027 +
                             -0.0001 * t0) * t0) * t0 +
                            (-0.0731 +
                             (0.0019 +
                              0.0009 * t0) * t0 +
                            (-0.0127 +
                              0.0011 * t0 + 0.0004 * t) * t) * t) * t) * t) * tas2r;

    // rotation matrix
    deuler("ZYZ", -zeta, theta, -z, mat);
}

}
