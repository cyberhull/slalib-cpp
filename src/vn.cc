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
 * Normalizes a 3-component vector, calculates vector's modulus (single precision).
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory. FORTRAN implementation was a procedure
 * that returned input vector modulus in one of its arguments; C++ implementation is a function returning the modulus
 * directly.
 *
 * @param vec A 3-component vector.
 * @param nvec Output: unit vector having the same direction as `vec`.
 *
 * @return Modulus of `vec`; if the modulus is zero, `nvec` is set to zero as well.
 */
float vn(const Vector<float> vec, Vector<float> nvec) {
    // calculate input vector's modulus
    float modulus = 0.0f;
    for (int i = 0; i < 3; i++) {
        float tmp = vec[i];
        modulus += tmp * tmp;
    }
    modulus = std::sqrt(modulus);
    float result = modulus;

    // normalize the vector
    if (modulus <= 0.0f) {
        modulus = 1.0f;
    }
    for (int i = 0; i < 3; i++) {
        nvec[i] = vec[i] / modulus;
    }
    return result;
}

}
