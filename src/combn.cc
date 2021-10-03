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
 * Generates the next combination, a subset of a specified size chosen from a specified number of items.
 *
 * This function returns, in the `list` array, a subset of `nsel` integers chosen from the range 1 to `ncand`
 * inclusive, in ascending order. Before calling the function for the first time, the caller must set the first
 * element of the `list` array to zero (any value less than 1 will do) to cause initialization.
 *
 * The first combination to be generated is:
 *
 *   list[0]=1, list[1]=2, ..., list[nsel-1]=nsel
 *
 * This is also the combination returned for the "finished" (`CPS_NO_MORE`) case. The final permutation to be
 * generated is:
 *
 *   list[0]=ncand, list[1]=ncand-1, ..., list[nsel-1]=ncand-nsel+1
 *
 * If the "finished" (`CPS_NO_MORE`) status is ignored, the function continues to deliver combinations, the pattern
 * repeating every ncand!/(nsel!*(ncand-nsel)!) calls.
 *
 * The algorithm is by R.F.Warren-Smith (private communication with P.T. Wallace).
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param nsel Number of items (subset size); must be at least 1; must be <= `ncand`.
 * @param ncand Number of candidates (set size); must be at least 1.
 * @param list Array of size `nsel`; `list[0..nsel-1]` represents latest combination; set `list[0]=0` to initialize.
 * @return Combination generation status, a `CPStatus` constant.
 */
CPStatus combn(int nsel, int ncand, int* list) {
    // validate arguments
    if (nsel < 1 || ncand < 1 || nsel > ncand) {
        return CPS_INVALID_ARG;
    }
    // just starting?
    if (list[0] < 1) {
        // yes: return 1,2,3...
        for (int j = 0; j < nsel;) {
            const int next_j = j + 1;
            list[j] = next_j;
            j = next_j;
        }
        return CPS_OK;
    }
    // no: find the first selection that we can increment
    // start with the first list item
    for (int i = 1;;) {
        // current list item.
        const int list_i = list[i - 1];

        // is this the final list item?
        const int nmax = (i >= nsel)?
            ncand + 1: // yes: comparison value is number of candidates plus one
            list[i]; // no: comparison value is next list item
        // can the current item be incremented?
        if (nmax - list_i > 1) {
            // yes: increment it.
            list[i - 1] = list_i + 1;

            // reinitialize the preceding items.
            for (int m = 0; m < i - 1;) {
                const int next_m = m + 1;
                list[m] = next_m;
                m = next_m;
            }
            return CPS_OK;
        } else {
            // can't increment the current item: is it the final one?
            if (i >= nsel) {
                // yes: restart the sequence
                for (int k = 0; k < nsel;) {
                    const int next_k = k + 1;
                    list[k] = next_k;
                    k = next_k;
                }
                return CPS_NO_MORE;
            } else {
                // no: proceed with next list item
                i++;
            }
        }
    }
}

}
