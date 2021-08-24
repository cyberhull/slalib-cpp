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
 * Calculates vector product of two 3-component vectors (double precision).
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param va First vector.
 * @param vb Second vector.
 * @param vc Output: vector product of the two first argument vectors; it is safe to specify either first of second
 *  vector as a result holder.
 */
void dvxv(const Vector<double> va, const Vector<double> vb, Vector<double> vc) {
    // form the vector product in a scratch vector
    // TODO: store first two expressions into scalars, store third into vc[2], copy scalars in vc[0] and vc[1]
    Vector<double> tmp;
    tmp[0]= va[1] * vb[2] - va[2] * vb[1];
    tmp[1]= va[2] * vb[0] - va[0] * vb[2];
    tmp[2]= va[0] * vb[1] - va[1] * vb[0];

    // return the result
    for (int i = 0; i < 3; i++) {
        vc[i] = tmp[i];
    }
}

}
