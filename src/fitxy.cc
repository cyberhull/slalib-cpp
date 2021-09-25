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
 * Fits a linear model to relate two sets of [x,y] coordinates.
 *
 * `sbr`, which is a boolean, selects the type of model fitted. Both `sbr` values produce a model `coeffs` which
 * consists of six coefficients, namely the zero points and, for each of XE and YE, the coefficient of XM and YM.
 * For `sbr`==`false`, all six coefficients are independent, modelling squash and shear as well as origin, scale,
 * and orientation. However, `sbr`==`true` selects the "solid body rotation" option; the model `coeffs` still
 * consists of the same six coefficients, but now two of them are used twice (appropriately signed). Origin, scale
 * and orientation are still modelled, but not squash or shear -- the units of X and Y have to be the same.
 *
 * For the returned `coeffs` (which can be fetched using `coeffs.get_a()` etc.), the model is:
 *
 *   XE = A + B*XM + C*YM
 *   YE = D + E*XM + F*YM
 *
 * For the "solid body rotation" option (`sbr`==`true`), the magnitudes of B and F, and of C and E, are equal. The
 * signs of these coefficients depend on whether there is a sign reversal between XE,YE and XM,YM;  fits are performed
 * with and without a sign reversal and the best one chosen.
 *
 * See also sla::pxy(), sla::invf(), sla::xy2xy(), and sla::dcmpf() functions.
 *
 * Original FORTRAN code by P.T. Wallace.
 *
 * @param sbr Whether to use "solid body rotation" (true), or 6 independent coefficients (false) model (boolean).
 * @param nsamples Number of samples; for "solid body rotation", `nsamples must be at least 2; for 6 independent
 *   coefficients, `nsamples` must be at least 3.
 * @param expected Expected [x,y] for each sample.
 * @param measured Measured [x,y] for each sample.
 * @param coeffs Return value: coefficients of the model.
 * @return A `FITResult` constant; if an error occurs and `FIT_INSUFFICIENT` is returned, `coeffs` are unchanged;
 *   if `FIT_NONE` is returned, then `coeffs` may have changed.
 */
FITStatus fitxy(bool sbr, int nsamples, const XYSamples expected, const XYSamples measured, FitCoeffs& coeffs) {
    int i, workspace[4];
    double x_expected, y_expected, x_measured, y_measured,
        sum_x_expected, sum_y_expected, sum_x_measured, sum_y_measured,
        vec[4], mat3[3][3], mat4[4][4], det, sdr2;

    // preset the status
    FITStatus status = FIT_OK;
    bool singular;

    // variable initializations to avoid compiler warnings
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
    double d = 0.0;
    double a_old = 0.0;
    double b_old = 0.0;
    double c_old = 0.0;
    double d_old = 0.0;
    double sdr2_old = 0.0;

    // float the number of samples
    const double num_samples = (double) nsamples;

    // 6 independent coefficients?
    if (sbr == false) {

        // six-coefficient linear model
        // ----------------------------

        // check enough samples
        if (nsamples >= 3) {

            // form summations
            sum_x_expected = 0.0;
            sum_y_expected = 0.0;
            sum_x_measured = 0.0;
            sum_y_measured = 0.0;
            double sum_xe_xm = 0.0;
            double sum_xe_ym = 0.0;
            double sum_ye_ym = 0.0;
            double sum_ye_xm = 0.0;
            double sum_xm_xm = 0.0;
            double sum_xm_ym = 0.0;
            double sum_ym_ym = 0.0;
            for (i = 0; i < nsamples; i++) {
                x_expected = expected[i][0];
                y_expected = expected[i][1];
                x_measured = measured[i][0];
                y_measured = measured[i][1];
                sum_x_expected += x_expected;
                sum_xe_xm += x_expected * x_measured;
                sum_xe_ym += x_expected * y_measured;
                sum_y_expected += y_expected;
                sum_ye_ym += y_expected * y_measured;
                sum_ye_xm += y_expected * x_measured;
                sum_x_measured += x_measured;
                sum_y_measured += y_measured;
                sum_xm_xm += x_measured * x_measured;
                sum_xm_ym += x_measured * y_measured;
                sum_ym_ym += y_measured * y_measured;
            }

            // solve for a,b,c in  x_expected = a + b*x_measured + c*y_measured
            vec[0] = sum_x_expected;
            vec[1] = sum_xe_xm;
            vec[2] = sum_xe_ym;
            mat3[0][0] = num_samples;
            mat3[0][1] = sum_x_measured;
            mat3[0][2] = sum_y_measured;
            mat3[1][0] = sum_x_measured;
            mat3[1][1] = sum_xm_xm;
            mat3[1][2] = sum_xm_ym;
            mat3[2][0] = sum_y_measured;
            mat3[2][1] = sum_xm_ym;
            mat3[2][2] = sum_ym_ym;
            singular = dmat(3, (double*) mat3, vec, det, workspace);
            if (singular == false) {
                for (i = 0; i < 3; i++) {
                    coeffs[i] = vec[i];
                }
                // solve for d,E,F in  y_expected = d + E*x_measured + F*y_measured
                vec[0] = sum_y_expected;
                vec[1] = sum_ye_xm;
                vec[2] = sum_ye_ym;
                dmxv(mat3, vec, coeffs.get_def());
            } else {
                // singular matrix, no 6-coefficient solution possible
                status = FIT_NONE;
            }
        } else {
            // insufficient data for 6-coefficient fit
            status = FIT_INSUFFICIENT;
        }
    } else {

        // four-coefficient solid body rotation model
        // ------------------------------------------

        // check enough samples
        if (nsamples >= 2) {

            // try two solutions, first without then with flip in X
            for (int solution = 1; solution <= 2; solution++) {
                const double sign = solution == 1? 1.0: -1.0;

                // form summations
                sum_x_expected = 0.0;
                sum_y_expected = 0.0;
                sum_x_measured = 0.0;
                sum_y_measured = 0.0;
                double sum_xx_yy = 0.0;
                double sum_xy_yx = 0.0;
                double sum_x2_y2 = 0.0;
                for (i = 0; i < nsamples; i++) {
                    x_expected = expected[i][0] * sign;
                    y_expected = expected[i][1];
                    x_measured = measured[i][0];
                    y_measured = measured[i][1];
                    sum_x_expected += x_expected;
                    sum_xx_yy += x_expected * x_measured + y_expected * y_measured;
                    sum_xy_yx += x_expected * y_measured - y_expected * x_measured;
                    sum_y_expected += y_expected;
                    sum_x_measured += x_measured;
                    sum_y_measured += y_measured;
                    sum_x2_y2 += x_measured * x_measured + y_measured * y_measured;
                }

                //
                // Solve for a,b,c,d in:  +/- x_expected = a + b*x_measured - c*y_measured
                //                          + y_expected = d + c*x_measured + b*y_measured
                vec[0] = sum_x_expected;
                vec[1] = sum_xx_yy;
                vec[2] = sum_xy_yx;
                vec[3] = sum_y_expected;
                mat4[0][0] = num_samples;
                mat4[0][1] = sum_x_measured;
                mat4[0][2] = -sum_y_measured;
                mat4[0][3] = 0.0;
                mat4[1][0] = sum_x_measured;
                mat4[1][1] = sum_x2_y2;
                mat4[1][2] = 0.0;
                mat4[1][3] = sum_y_measured;
                mat4[2][0] = sum_y_measured;
                mat4[2][1] = 0.0;
                mat4[2][2] = -sum_x2_y2;
                mat4[2][3] = -sum_x_measured;
                mat4[3][0] = 0.0;
                mat4[3][1] = sum_y_measured;
                mat4[3][2] = sum_x_measured;
                mat4[3][3] = num_samples;
                singular = dmat(4, (double*) mat4, vec, det, workspace);
                if (singular == false) {
                    a = vec[0];
                    b = vec[1];
                    c = vec[2];
                    d = vec[3];

                    // determine sum of radial errors squared
                    sdr2 = 0.0;
                    for (i = 0; i < nsamples; i++) {
                        x_measured = measured[i][0];
                        y_measured = measured[i][1];
                        const double xr = a + b * x_measured - c * y_measured - expected[i][0] * sign;
                        const double yr = d + c * x_measured + b * y_measured - expected[i][1];
                        sdr2 = sdr2 + xr * xr + yr * yr;
                    }
                } else {
                    // singular: set flag
                    sdr2 = -1.0;
                }

                // if first pass and non-singular, save variables
                if (solution == 1 && singular == false) {
                    a_old = a;
                    b_old = b;
                    c_old = c;
                    d_old = d;
                    sdr2_old = sdr2;
                }

            }

            // pick the best of the two solutions
            if (sdr2_old >= 0.0 && (sdr2_old <= sdr2 || nsamples == 2)) {
                coeffs[0] = a_old;
                coeffs[1] = b_old;
                coeffs[2] = -c_old;
                coeffs[3] = d_old;
                coeffs[4] = c_old;
                coeffs[5] = b_old;
            } else if (singular == false) {
                coeffs[0] = -a;
                coeffs[1] = -b;
                coeffs[2] = c;
                coeffs[3] = d;
                coeffs[4] = c;
                coeffs[5] = b;
            } else {
                // no 4-coefficient fit possible
                status = FIT_NONE;
            }
        } else {

            // insufficient data for 4-coefficient fit
            status = FIT_INSUFFICIENT;
        }
    }
    return status;
}

}
