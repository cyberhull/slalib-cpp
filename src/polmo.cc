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
 * Corrects site longitude and latitude for polar motion and calculates azimuth difference between celestial and
 * terrestrial poles (double precision).
 *
 * "Mean" longitude and latitude are the (fixed) values for the site's location with respect to the IERS terrestrial
 * reference frame; the latitude is geodetic. TAKE CARE WITH THE LONGITUDE SIGN CONVENTION. The longitudes used by
 * the present function are east-positive, in accordance with geographical convention (and right-handed). In
 * particular, note that the longitudes returned by the sla::obs() function are west-positive, following astronomical
 * usage, and must be reversed in sign before use in the present function.
 *
 * `x_pm` and `y_pm` are the (changing) coordinates of the Celestial Ephemeris Pole with respect to the IERS
 * Reference Pole. `x_pm` is positive along the meridian at longitude 0 degrees, and `y_pm` is positive along the
 * meridian at longitude 270 degrees (i.e. 90 degrees west). Values for `x_pm`, `y_pm` can be obtained from IERS
 * circulars and equivalent publications; the maximum amplitude observed so far is about 0.3 arcseconds.
 *
 * "True" longitude and latitude are the (moving) values for the site's location with respect to the celestial
 * ephemeris pole and the meridian which corresponds to the Greenwich apparent sidereal time. The true longitude and
 * latitude link the terrestrial coordinates with the standard celestial models (for precession, nutation, sidereal
 * time etc).
 *
 * The azimuths produced by sla::aop() and sla::aopqk() are with respect to due north as defined by the Celestial
 * Ephemeris Pole, and can therefore be called "celestial azimuths". However, a telescope fixed to the Earth measures
 * azimuth essentially with respect to due north as defined by the IERS Reference Pole, and can therefore be called
 * "terrestrial azimuth". Uncorrected, this would manifest itself as a changing "azimuth zero-point error". The
 * value `d_az` is the correction to be added to a celestial azimuth to produce a terrestrial azimuth.
 *
 * The present function is rigorous. For most practical purposes, the following simplified formulae provide an
 * adequate approximation:
 *
 *   `t_long` = `m_long` + `x_pm` * cos(`m_long`) - `y_pm` * sin(`m_long`)
 *   `t_phi`  = `m_phi` + (`x_pm` * sin(`m_long`) + `y_pm` * cos(`m_long`)) * tan(`m_phi`)
 *   `d_az`   = -sqrt(`x_pm` * `x_pm` + `y_pm` * `y_pm`) * cos(`m_long` - atan2(`x_pm`,`y_pm`)) / cos(`m_phi`)
 *
 * An alternative formulation for `d_az` is:
 *
 *   x = cos(`m_long`) * cos(`m_phi`)
 *   y = sin(`m_long`) * cos(`m_phi`)
 *   `d_az` = atan2(-x * `y_pm` - y * `x_pm`, x * x + y * y)
 *
 * Reference:
 *   Seidelmann, P.K. (ed), 1992. "Explanatory Supplement to the Astronomical Almanac", ISBN 0-935702-68-7,
 *   sections 3.27, 4.25, 4.52.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param m_long Mean longitude of the observer (radians, east positive).
 * @param m_phi Mean geodetic latitude of the observer (radians).
 * @param x_pm Polar motion x-coordinate (radians).
 * @param y_pm Polar motion y-coordinate (radians).
 * @param t_long Return value: true longitude of the observer (radians, east positive).
 * @param t_phi Return value: true geodetic latitude of the observer (radians).
 * @param d_az Return value: azimuth correction (terrestrial-celestial, radians).
 */
void polmo(double m_long, double m_phi, double x_pm, double y_pm, double& t_long, double& t_phi, double& d_az) {
    // site mean longitude and mean geodetic latitude as a Cartesian vector
    double sin_long = std::sin(m_long);
    double cos_long = std::cos(m_long);
    const double sin_phi = std::sin(m_phi);
    double cos_phi = std::cos(m_phi);

    const double xm = cos_long * cos_phi;
    const double ym = sin_long * cos_phi;
    const double zm = sin_phi;

    // rotate site vector by polar motion, Y-component then X-component
    const double sin_xpm = std::sin(x_pm);
    const double cos_xpm = std::cos(x_pm);
    const double sin_ypm = std::sin(y_pm);
    const double cos_ypm = std::cos(y_pm);

    const double zw = (-ym * sin_ypm + zm * cos_ypm);

    double xt = xm * cos_xpm - zw * sin_xpm;
    const double yt = ym * cos_ypm + zm * sin_ypm;
    const double zt = xm * sin_xpm + zw * cos_xpm;

    // rotate also the geocentric direction of the terrestrial pole (0,0,1)
    const double xnm = -sin_xpm * cos_ypm;
    const double ynm = sin_ypm;
    const double znm = cos_xpm * cos_ypm;

    cos_phi = std::sqrt(xt * xt + yt * yt);
    if (cos_phi == 0.0) {
        xt = 1.0;
    }
    sin_long = yt / cos_phi;
    cos_long = xt / cos_phi;

    // return true longitude and true geodetic latitude of site
    t_long = (xt != 0.0 || yt != 0.0)? std::atan2(yt, xt): 0.0;
    t_phi = std::atan2(zt, cos_phi);

    // return current azimuth of terrestrial pole seen from site position
    const double xnt = (xnm * cos_long + ynm * sin_long) * zt - znm * cos_phi;
    const double ynt = -xnm * sin_long + ynm * cos_long;
    d_az = (xnt != 0.0 || ynt != 0.0)? std::atan2(-ynt, -xnt): 0.0;
}

}
