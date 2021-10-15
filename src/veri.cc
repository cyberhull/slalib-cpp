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
 * Reports the SLALIB version number as an integer.
 *
 * To obtain the version number in a printable form, see function sla::vers().
 *
 * Version numbers for the C++ SLALIB port had been taken from the `component id="sla"` tag of the `componentset.xml`
 * file located in the root directory of the `starlink` project.
 *
 * Original FORTRAN code by Norman Gray / Council for the Central Laboratory of the Research Councils.
 *
 * @return SLALIB version number as an integer `M` * 1,000,000 + `N` * 1,000 + `R`, where `M` is the major version,
 *   `N` the minor version, and `R` the release number.
 */
int veri() {
    return SLALIB_PACKAGE_VERSION_MAJOR * 1000000 + SLALIB_PACKAGE_VERSION_MINOR * 1000 +
        SLALIB_PACKAGE_VERSION_RELEASE;
}

}
