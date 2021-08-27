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
  av2m.cc dav2m.cc cc2s.cc dcc2s.cc cs2c.cc dcs2c.cc euler.cc deuler.cc imxv.cc dimxv.cc
  m2av.cc dm2av.cc mxm.cc dmxm.cc mxv.cc dmxv.cc vdv.cc dvdv.cc vn.cc dvn.cc vxv.cc dvxv.cc
  zd.cc pa.cc
  bear.cc dbear.cc pav.cc dpav.cc
  e2h.cc de2h.cc h2e.cc dh2e.cc
  caf2r.cc daf2r.cc
  cldj.cc caldj.cc clyd.cc calyd.cc djcal.cc djcl.cc
  cd2tf.cc dd2tf.cc cr2af.cc dr2af.cc cr2tf.cc dr2tf.cc
  ctf2d.cc dtf2d.cc ctf2r.cc dtf2r.cc
  dat.cc dt.cc dtt.cc rcc.cc
  range.cc drange.cc ranorm.cc dranrm.cc
  atmdsp.cc refcoq.cc refro.cc refco.cc refv.cc refz.cc
  ecmat.cc dmat.cc smat.cc altaz.cc
  nutc.cc nut.cc nutc80.cc
  epj2d.cc epj.cc epb2d.cc epb.cc
  prec.cc precl.cc prenut.cc
  dsepv.cc sepv.cc dsep.cc sep.cc
  prebn.cc preces.cc supgal.cc
  rverot.cc rvgalc.cc rvlg.cc rvlsrd.cc rvlsrk.cc
  cc62s.cc dc62s.cc
  wait.cc
  f77_utils.h
  slalib.cc slalib.h)
