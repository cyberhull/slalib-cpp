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
 * Given arrays of "expected" and "measured" [x,y] coordinates, and a linear model relating them (as produced by
 * sla::fitxy()), computes the array of "predicted" coordinates and the RMS residuals (double precision).
 *
 * The model is supplied in the array/structure `model` (with coefficients accessible as `model.get_a()`,
 * `model.get_b(), etc.); the model is applied as follows (suffix 'P' stands for "predicted", 'M' for "measured",
 * and 'E' for "expected"):
 *
 *   XP = A + B * XM + C * YM
 *   YP = D + E * XM + F * YM
 *
 * The residuals are (XP-XE) and (YP-YE).
 *
 * See also sla::fitxy(), sla::invf(), sla::xy2xy(), and sla::dcmpf() functions.
 *
 * Original FORTRAN code by P.T. Wallace / Rutherford Appleton Laboratory.
 *
 * @param nsamples Number of samples; If `nsamples` is less than or equal to zero, no coordinates are transformed,
 *   and the RMS residuals are all zero.
 * @param expected Expected [x,y] for each sample.
 * @param measured Measured [x,y] for each sample.
 * @param model Coefficients of model.
 * @param predicted Predicted [x,y] for each sample.
 * @param x_rms Return value: RMS in X.
 * @param y_rms Return value: RMS in Y.
 * @param rms Return value: total RMS (vector sum of `x_rms` and `y_rms`).
 */
void pxy(int nsamples, const XYSamples expected, const XYSamples measured, const FitCoeffs& model,
    XYSamples predicted, double& x_rms, double& y_rms, double& rms) {

    // initialize summations
    double sum_dx2 = 0.0;
    double sum_dy2 = 0.0;

    // loop by sample
    for (int i = 0; i < nsamples; i++) {
        // transform "measured" [X,Y] to "predicted" [X,Y]
        double x_predicted, y_predicted;
        xy2xy(measured[i][0], measured[i][1], model, x_predicted, y_predicted);
        predicted[i][0] = x_predicted;
        predicted[i][1] = y_predicted;

        // compute residuals in X and Y, and update summations
        const double dx = expected[i][0] - x_predicted;
        const double dy = expected[i][1] - y_predicted;
        sum_dx2 += dx * dx;
        sum_dy2 += dy * dy;
    }

    // compute RMS values
    const double p = std::max(1.0, (double) nsamples);
    x_rms = std::sqrt(sum_dx2 / p);
    y_rms = std::sqrt(sum_dy2 / p);
    rms = std::sqrt(x_rms * x_rms + y_rms * y_rms);
}

}
