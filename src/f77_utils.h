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
#ifndef SLALIB_F77_UTILS_H_INCLUDED
#define SLALIB_F77_UTILS_H_INCLUDED

#include <cmath>
#include <type_traits>

namespace sla {

// C++ counterparts of the F77 math functions; to be inlined at a later stage
template <typename T, std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
inline T f_sign(T a, T b) {
    const T aa = std::abs(a);
    return b >= 0.0? aa: -aa;
}

template <typename T, std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
inline int f_nint(T n) { return (int) std::round(n); }

template <typename T, std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
inline double f_aint(T n) { return std::trunc(n); }

template <typename T, std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
inline double f_anint(T n) { return std::round(n); }

} // sla

#endif // SLALIB_F77_UTILS_H_INCLUDED
