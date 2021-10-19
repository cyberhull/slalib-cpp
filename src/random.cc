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
#include "f77_utils.h"
#include <cmath>
#include <random>

namespace sla {

/**
 * Generates pseudo-random real number in the range [0.0f..1.0f) (single precision).
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param seed An arbitrary `float` number; used first time through *only*.
 * @return Pseudo-random `float` number in the range [0.0f..1.0f).
 */
float random(float seed) {
    static SLALIB_THREAD_LOCAL std::default_random_engine generator;
    static SLALIB_THREAD_LOCAL bool first_time = true;
    if (first_time) {
        /*
         * The following numeric gymnastics is, strictly speaking, not needed: we could have used some much simpler
         * conversion of the floating-point `seed` to an integer one (original FORTRAN implementation took care
         * not to pass `1` seed to `RAND()` [by converting `seed` to a "large integer"] as that would mean
         * "initialize / restart" the generator). It's been decided to port the original implementation for maximum
         * compatibility of seeding the generator.
         */
        float aseed = std::fabs(seed) + 1.0f;
        int iseed = f_nint(aseed / std::pow(10.0f, (f_nint(std::log10(aseed)) - 6)));
        if ((iseed & 1) == 0) {
            iseed++;
        }
        generator.seed(iseed);
        first_time = false;
    }
    return std::generate_canonical<float, std::numeric_limits<float>::digits>(generator);
}

}
