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
 * Calculates angle between two vectors (single precision).
 *
 * The simplest formulation would use dot product alone. However, this would reduce the accuracy for angles near zero
 * and pi. The algorithm uses both cross product and dot product, which maintains accuracy for all sizes of angle.
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param va First vector; does not have to be unit length; if null, zero is returned.
 * @param vb Second vector; does not have to be unit length; if null, zero is returned.
 * @return The angle, in radians, between the two vectors; it is always positive.
 */
float sepv(const Vector<float> va, const Vector<float> vb) {
    // use double precision version.
    const Vector<double> dv1 = {va[0], va[1], va[2]};
    const Vector<double> dv2 = {vb[0], vb[1], vb[2]};
    return (float) dsepv(dv1, dv2);
}

}
