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
 * Normalizes a 3-component vector, calculates vector's modulus (double precision).
 *
 * Original FORTRAN code by P.T. Wallace. FORTRAN implementation was a procedure that returned input vector modulus
 * in one of its arguments; C++ implementation is a function returning the modulus directly.
 *
 * @param v A 3-component vector.
 * @param uv Output: unit vector having the same direction as `v`.
 *
 * @return Modulus of `v`; if the modulus is zero, `uv` is set to zero as well.
 */
double dvn(const vector<double> v, vector<double> uv) {
    // calculate input vector's modulus
    double modulus = 0.0;
    for (int i = 0; i < 3; i++) {
        double tmp = v[i];
        modulus += tmp * tmp;
    }
    modulus = std::sqrt(modulus);
    double result = modulus;

    // normalize the vector
    if (modulus <= 0.0) {
        modulus = 1.0;
    }
    for (int i = 0; i < 3; i++) {
        uv[i] = v[i] / modulus;
    }
    return result;
}

}
