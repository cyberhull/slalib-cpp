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
 * @param v1 First vector; does not have to be unit length; if null, zero is returned.
 * @param v2 Second vector; does not have to be unit length; if null, zero is returned.
 * @return The angle, in radians, between the two vectors; it is always positive.
 */
float sepv(const Vector<float> v1, const Vector<float> v2) {
    // use double precision version.
    const Vector<double> dv1 = {v1[0], v1[1], v1[2]};
    const Vector<double> dv2 = {v2[0], v2[1], v2[2]};
    return (float) dsepv(dv1, dv2);
}

}
