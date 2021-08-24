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
 * Converts Modified Julian Date to Julian Epoch (double precision).
 *
 * Reference:
 *   Lieske,J.H., 1979. Astron.Astrophys., 73, 282.
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param mjd Modified Julian Date (JD - 2400000.5).
 * @return Julian Epoch.
 */
double epj(double mjd) {
    return 2000.0 + (mjd - 51544.5) / 365.25;
}

}
