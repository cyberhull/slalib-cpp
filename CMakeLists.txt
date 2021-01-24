#
# C++ Port of the SLALIB library.
# Written by Vadim Sytnikov.
# Copyright (C) 2021 CyberHULL, Ltd.
# All rights reserved.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# -----------------------------------------------------------------------------
#
# Main project file.
#
cmake_minimum_required(VERSION 3.10)
project(SLALibCpp)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_C_OUTPUT_EXTENSION ".o")
set(CMAKE_CXX_OUTPUT_EXTENSION ".o")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# use "-s -O3" for optimized builds
set(SLA_COMMON_FLAGS "-Wall")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SLA_COMMON_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SLA_COMMON_FLAGS} -std=c++14")

include_directories(src)

add_subdirectory(src)
add_subdirectory(tests)
