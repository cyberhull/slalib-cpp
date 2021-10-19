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
#include <random>

namespace sla {

/**
 * Generates pseudo-random normal deviate ( = 'Gaussian residual') (single precision).
 *
 * The Box-Muller algorithm is used. This is described in Numerical Recipes, section 7.2.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param stdev Standard deviation.
 * @return The results of many calls to this function will be normally distributed with mean zero and standard
 *   deviation `stdev`.
 */
float gresid(float stdev) {
    static SLALIB_THREAD_LOCAL std::default_random_engine generator;
    static SLALIB_THREAD_LOCAL float g_next;
    static SLALIB_THREAD_LOCAL bool ftf = true, first = true;

    auto next_rn = []() -> float {
        return std::generate_canonical<float, std::numeric_limits<float>::digits>(generator);
    };

    float g;

    // first time through, initialise the random-number generator
    if (ftf) {
        generator.seed(123456789);
        ftf = false;
    }

    // second normal deviate of the pair available?
    if (first) {

        // no - generate two random numbers inside unit circle
        float x, y, r;
        do {
            // generate two random numbers in range +/- 1
            x = 2.0f * next_rn() - 1.0f;
            y = 2.0f * next_rn() - 1.0f;

            // try again if not in unit circle
            r = x * x + y * y;
        } while (r >= 1.0f);

        // Box-Muller transformation, generating two deviates
        const float w = std::sqrt(-2.0f * std::log(r) / std::max(r, 1.0e-20f));
        g_next = x * w;
        g = y * w;

        // set flag to indicate availability of next deviate
        first = false;
    } else {
        // return second deviate of the pair & reset flag
        g = g_next;
        first = true;
    }

    // scale the deviate by the required standard deviation
    return g * stdev;
}

}
