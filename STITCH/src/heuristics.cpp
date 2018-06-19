// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include <iomanip>
#include <string>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
using namespace Rcpp;



//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_make_smoothed_rate(
    const Rcpp::NumericVector & sigmaSum_unnormalized,
    const Rcpp::NumericVector & sigma_rate,
    const Rcpp::IntegerVector & L_grid,
    const Rcpp::IntegerVector & grid_distances,
    const int nGrids,
    const int shuffle_bin_radius
) {

    int iSNP, iSNP_left, iSNP_right, bp_remaining, bp_to_add;
    int focal_point;
    int min_L_grid = Rcpp::min(L_grid);
    int max_L_grid = Rcpp::max(L_grid);
    Rcpp::NumericVector smoothed_rate(nGrids - 1);
    for(iSNP = 0; iSNP < (nGrids - 1); iSNP++) {
        focal_point = ((L_grid(iSNP) + L_grid(iSNP + 1)) / 2); // automatically rounded
        if (
            (min_L_grid <= (focal_point - shuffle_bin_radius)) &
            ((focal_point + shuffle_bin_radius) <= max_L_grid)
            ) {
            // left
            iSNP_left=iSNP;
            bp_remaining = shuffle_bin_radius;
            while (0 < bp_remaining) {
                bp_to_add = (L_grid(iSNP_left + 1) - L_grid(iSNP_left) + 1);
                if ((bp_remaining - bp_to_add) < 0) {
                    bp_to_add = bp_remaining;
                    bp_remaining = 0;
                } else {
                    bp_remaining = bp_remaining - bp_to_add;
                    iSNP_left = iSNP_left - 1;
                }
                // add bit
                smoothed_rate(iSNP) = smoothed_rate(iSNP) +     \
                    bp_to_add * sigma_rate(iSNP_left);
            }
            // right
            iSNP_right = iSNP + 1;
            bp_remaining = shuffle_bin_radius;
            while(0 < bp_remaining) {
                bp_to_add = (L_grid(iSNP_right) - L_grid(iSNP_right - 1) + 1);
                if ((bp_remaining - bp_to_add) < 0) {
                    bp_to_add = bp_remaining;
                    bp_remaining = 0;
                } else {
                    bp_remaining = bp_remaining - bp_to_add;
                    iSNP_right = iSNP_right + 1;
                }
                // add bit
                smoothed_rate(iSNP) = smoothed_rate(iSNP) +     \
                    bp_to_add * sigma_rate(iSNP_right);
            }
        } else {
            smoothed_rate(iSNP) = NA_REAL;
        }
    }
    
    return(smoothed_rate);
    
}
