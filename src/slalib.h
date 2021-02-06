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
#ifndef _SLALIB_H
#define _SLALIB_H

namespace sla {

template<typename T>
using vector = T[3];

template<typename T>
using matrix = T[3][3];

double airmas(double zenith_dist);
void cc2s(const vector<float> vec, matrix<float> mat);

} // sla

#endif // _SLALIB_H
