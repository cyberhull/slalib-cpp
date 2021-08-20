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
 * Applies precession - either FK4 (Bessel-Newcomb, pre IAU 1976) or FK5 (Fricke, post IAU 1976) as required.
 *
 * The epochs are Besselian if `system`==`CAT_FK4` and Julian if it's `CAT_FK5`. For example, to precess coordinates
 * in the old system from equinox 1900.0 to 1950.0 the call would be:
 *   SphericalDir<double> pos = {<RA>, <Dec>};
 *   preces(CAT_FK4, 1900.0, 1950.0, pos);
 *
 * This function will NOT correctly convert between the old and the new systems - for example conversion from B1950
 * to J2000. For these purposes see sla::fk425(), sla::fk524(), sla::fk45z(), and sla::fk54z().
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param system Precession to be applied, a `CAT_xxx` constant (only `CAT_FK4` or `CAT_FK5` are accepted).
 * @param ep0 Starting epoch.
 * @param ep1 Ending epoch.
 * @param pos RA,Dec: mean equator & equinox of epoch `ep0` (argument) or epoch `ep1` (return value); if an invalid
 *   `system` is specified, values of -99.0,-99.0 will be returned.
 */
void preces(Catalogue system, double ep0, double ep1, SphericalDir<double>& pos) {
    // validate reference system
    switch (system) {
        case CAT_FK4:
        case CAT_FK5:
            // generate appropriate precession matrix
            matrix<double> m_precession;
            if (system == CAT_FK4) {
                prebn(ep0, ep1, m_precession);
            } else {
                prec(ep0, ep1, m_precession);
            }
            // convert ra,Dec to x,y,z
            vector<double> v1;
            dcs2c(pos, v1);

            // precess
            vector<double> v2;
            dmxv(m_precession, v1, v2);

            // convert back to ra,Dec
            dcc2s(v2, pos);
            pos.sd_a = dranrm(pos.sd_a);
            break;
        default:
            pos.sd_a = -99.0;
            pos.sd_b = -99.0;
    }
}

}
