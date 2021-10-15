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

/*
 * The macro to be used to stringify non-string values is `SLALIB_STRINGIFY()`; the `SLALIB_STRINGIFY_HELPER()` macro
 * is there to supply correct string to the `#` operator; without it, `SLALIB_STRINGIFY(name) #name` would always
 * return string "name" (i.e. formal argument in quotes).
 */
#define SLALIB_STRINGIFY(name) SLALIB_STRINGIFY_HELPER(name)
#define SLALIB_STRINGIFY_HELPER(str) #str
#define SLALIB_VERSION_STRING SLALIB_STRINGIFY(SLALIB_PACKAGE_VERSION_MAJOR) "." \
    SLALIB_STRINGIFY(SLALIB_PACKAGE_VERSION_MINOR) "-" SLALIB_STRINGIFY(SLALIB_PACKAGE_VERSION_RELEASE)

/**
 * Reports the SLALIB version number as a string.
 *
 * To obtain the version number in a more easily processed form, see function sla::veri().
 *
 * Version numbers for the C++ SLALIB port had been taken from the `component id="sla"` tag of the `componentset.xml`
 * file located in the root directory of the `starlink` project.
 *
 * Original FORTRAN code by Norman Gray / Council for the Central Laboratory of the Research Councils.
 *
 * @return SLALIB version number, in the form "M.N-R"; the major version is `M`, the minor version is `N`, and
 *   release is `R`.
 */
const char* vers() {
    return SLALIB_VERSION_STRING;
}

}
