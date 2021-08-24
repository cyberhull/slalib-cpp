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
 * Calculates angle between two vectors (double precision).
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
double dsepv(const Vector<double> v1, const Vector<double> v2) {
    // modulus of cross product = sine multiplied by the two moduli
    Vector<double> v1_x_v2, wv;
    dvxv(v1, v2, v1_x_v2);
    const double s = dvn(v1_x_v2, wv);

    // dot product = cosine multiplied by the two moduli.
    const double c = dvdv(v1, v2);

    // angle between the vectors
    return s != 0.0 || c != 0.0? std::atan2(s,c): 0.0;
}

}
