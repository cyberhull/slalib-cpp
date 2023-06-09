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
    dat.cc dt.cc dtt.cc rcc.cc gmst.cc gmsta.cc
    range.cc drange.cc ranorm.cc dranrm.cc
    atmdsp.cc refcoq.cc refro.cc refco.cc refv.cc refz.cc
    ecmat.cc dmat.cc smat.cc svd.cc svdsol.cc svdcov.cc altaz.cc
    nutc.cc nut.cc nutc80.cc
    epj2d.cc epj.cc epb2d.cc epb.cc epco.cc
    prec.cc precl.cc prenut.cc
    dsepv.cc sepv.cc dsep.cc sep.cc
    prebn.cc preces.cc supgal.cc
    rverot.cc rvgalc.cc rvlg.cc rvlsrd.cc rvlsrk.cc
    cc62s.cc dc62s.cc cs2c6.cc ds2c6.cc
    etrms.cc addet.cc subet.cc
    geoc.cc pvobs.cc pcd.cc unpcd.cc
    eqeqx.cc eqecl.cc eqgal.cc galeq.cc
    fitxy.cc xy2xy.cc pxy.cc invf.cc dcmpf.cc
    pm.cc earth.cc ecor.cc ecleq.cc polmo.cc galsup.cc
    s2tp.cc ds2tp.cc tp2s.cc dtp2s.cc tps2c.cc dtps2c.cc
    tp2v.cc dtp2v.cc v2tp.cc dv2tp.cc tpv2c.cc dtpv2c.cc
    combn.cc permut.cc
    evp.cc epv.cc
    eg50.cc ge50.cc
    pdq2h.cc pda2h.cc
    veri.cc vers.cc random.cc gresid.cc wait.cc
    moon.cc dmoon.cc
    obs.cc
    f77_utils.h
    slalib.cc slalib.h)
