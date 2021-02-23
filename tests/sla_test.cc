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
#include <cstdio>
#include <cstring>
#include <cmath>
#include <resolv.h>

#include "sla_test.h"
#include "../src/slalib.h"

namespace sla {

///////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS
///////////////////////////////////////////////////////////////////////////////

// reports error in a function test
static void err(const char* func, const char* test, bool &status) {
    std::printf("Test '%s' of the function %s() FAILED\n", test, func);
    status = false;
}

/// Validates a character string result.
static void vcs(const char* str, const char* str_ok, const char* func, const char* test, bool &status) {
    if (std::strcmp(str, str_ok) != 0) {
        err(func, test, status);
        std::printf("\tExpected: '%s'\n", str_ok);
        std::printf("\t  Actual: '%s'\n", str);
    }
}

/// Validates an integer result.
static void viv(const int val, const int val_ok, const char* func, const char* test, bool &status) {
    if (val != val_ok) {
        err(func, test, status);
        std::printf("\tExpected: %d\n", val_ok);
        std::printf("\t  Actual: %d\n", val);
    }
}

/// Validates a long result; in FORTRAN implementation of SLALIB, "long" values are INTEGER*4, so C/C++ `int`s
/// work just fine, there is no need to use `long`s.
static void vlv(const int val, const int val_ok, const char* func, const char* test, bool &status) {
    static_assert(sizeof val >= 4, "`int` values must be at least 32-bit");
    viv(val, val_ok, func, test, status);
}

/// Validates a double precision floating point result.
static void vvd(const double val, const double val_ok, const double tolerance,
    const char* func, const char* test, bool &status) {
    if (std::fabs(val - val_ok) > tolerance)  {
        err(func, test, status);
        std::printf("\tExpected: %30.19f\n", val_ok);
        std::printf("\t  Actual: %30.19f\n", val);
    }
}

///////////////////////////////////////////////////////////////////////////////
// INDIVIDUAL FUNCTION TESTS
///////////////////////////////////////////////////////////////////////////////

/// Tests sla::airmas() function.
void t_airmas(bool& status) {
    vvd(airmas(1.2354), 3.015698990074724, 1e-12, "sla::airmas", "", status);
}

///////////////////////////////////////////////////////////////////////////////
// MODULE ENTRY POINT
///////////////////////////////////////////////////////////////////////////////

/// Tests all SLALIB functions and procedures.
bool sla_test() {
    bool status = true;
    t_airmas(status);
    return status;
}

}
