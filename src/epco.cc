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
 * Converts an epoch into the appropriate form - Besselian or Julian (double precision).
 *
 * The result is always either equal to or very close to the given epoch `epoch`. The function is required only in
 * applications where punctilious treatment of heterogeneous mixtures of star positions is necessary.
 *
 * @param result Type of the result: 'B', 'b' (Besselian), 'J', or 'j' (Julian).
 * @param given Type of the given epoch: 'B', 'b' (Besselian), 'J', or 'j' (Julian).
 * @param epoch Epoch to convert.
 * @return If `result` and `given` are the same, returns specified `epoch`; otherwise, if either `result` or `given`
 *   is invalid (not 'B', 'b', 'J', or 'j'), returns 0.0; otherwise, returns converted (B to J, or J to B) value.
 */
double epco(char result, char given, double epoch) {
    if (result == 'b') {
        result = 'B';
    }
    if (result != 'B' && result != 'J') {
        return 0.0;
    }
    if (given == 'j') {
        given = 'J';
    }
    if (given != 'B' && given != 'J') {
        return 0.0;
    }
    if (result == given) {
        return epoch;
    } else if (result == 'B') {
        return epb(epj2d(epoch));
    } else { // 'J'
        return epj(epb2d(epoch));
    }
}

}
