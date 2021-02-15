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
# The library modules.
#
project(SLALibrary)

add_library(slalib
  airmas.cc
  av2m.cc dav2m.cc
  cc2s.cc dcc2s.cc
  cs2c.cc dcs2c.cc
  euler.cc deuler.cc
  imxv.cc dimxv.cc
  m2av.cc dm2av.cc
  mxm.cc
  slalib.cc slalib.h)
