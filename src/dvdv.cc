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
 * Computes scalar product of two 3-element vectors (double precision).
 *
 * @param va First vector.
 * @param vb Second vector.
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @return Scalar product va.vb
 */
double dvdv(const vector<double> va, const vector<double> vb) {
    return va[0] * vb[0] + va[1] * vb[1] + va[2] * vb[2];
}

}
