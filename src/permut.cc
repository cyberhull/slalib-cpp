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
 * Generates the next permutation of a specified number of items.
 *
 * This function returns, in the `order[]` array, the integers 1 to N inclusive, in an order that depends on the
 * current contents of the `state[]` array. Before calling the function for the first time, the caller must set the
 * first (0 index) element of the `state[]` array to -1 (any negative number will do) to cause the `state[]` array
 * to be fully initialized.
 *
 * The first permutation to be generated is:
 *
 *   order[0]=n, order[1]=n-1, ..., order[n-1]=1
 *
 * This is also the permutation returned for the "finished" (`CPS_NO_MORE`) case. The final permutation to be
 * generated is:
 *
 *   order[0]=1, iorder[1]=2, ..., iorder[n-1]=n
 *
 * If the "finished" (`CPS_NO_MORE`) status is ignored, the routine continues to deliver permutations, the pattern
 * repeating every n! calls.
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param n Number of items: there will be N! permutations.
 * @param state Array of size `n`, the state; set `state[0]=-1` to initialize.
 * @param order Next permutation of numbers 1,2,...,N.
 * @return Permutation generation status, a `CPStatus` constant.
 */
CPStatus permut(int n, int* state, int* order) {
    // validate arguments
    if (n < 1) {
        return CPS_INVALID_ARG;
    }

    // if just starting, initialize state array
    if (state[0] < 0) {
        state[0] = -1;
        for (int i = 1; i < n; i++) {
            state[i] = 0;
      }
    }

    /*
     *  Increment the state number
     *  ---------------------------------------------------------------------------
     *  The state number, maintained in the state array, is a mixed-radix number
     *  with n! states. The least significant digit, with a radix of 1, is in
     *  state[0]. The next digit, in state[1], has a radix of 2, and so on.
     */
    // increment the least-significant digit of the state number
    state[0]++;

    // digit by digit starting with the least significant
    CPStatus status = CPS_OK;
    for (int j = 1; j <= n; j++) {
        // carry?
        if (state[j - 1] >= j) {
            // yes: reset the current digit
            state[j - 1] = 0;
            // overflow?
            if (j >= n) {
                // yes: there are no more permutations
                status = CPS_NO_MORE;
            } else {
                // no: carry
                state[j]++;
            }
        }
    }

    /*
     * Translate the state number into the corresponding permutation order
     * ------------------------------------------------------------------------
     */

    // initialize the order array; all but one element will be overwritten
    for (int k = 0; k < n; k++) {
        order[k] = 1;
    }

    // look at each state number digit, starting with the most significant
    for (int l = n; l >= 2; l--) {
        // initialize the position where the new number will go
        int slot = 0;
        // the state number digit says which unfilled slot is to be used
        for (int skip = 0; skip <= state[l - 1]; skip++) {
            // increment the slot number until an unused slot is found
            slot++;
            while (order[slot - 1] > 1) {
                slot++;
            }
        }
        // store the number in the permutation order array.
        order[slot - 1] = l;
    }
    return status;
}

}
