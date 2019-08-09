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
    const Rcpp::NumericVector & sigma_rate,
    const Rcpp::IntegerVector & L_grid,
    const int shuffle_bin_radius,
    const bool verbose = false
) {
    const int nGrids = L_grid.length();
    int iGrid, iGrid_left, iGrid_right, bp_remaining, bp_to_add;
    int focal_point, bp_prev;
    double total_bp_added;
    //int min_L_grid = Rcpp::min(L_grid);
    //int max_L_grid = Rcpp::max(L_grid);
    Rcpp::NumericVector smoothed_rate(nGrids - 1);
    for(iGrid = 0; iGrid < (nGrids - 1); iGrid++) {
        if (verbose) {
            std::cout << "iGrid=" << iGrid << std::endl;
        }
        focal_point = (L_grid(iGrid) + L_grid(iGrid + 1)) / 2; // automatically rounded        
        if (verbose) {
            std::cout << "focal_point=" << focal_point << std::endl;
        }
        //
        // left
        //
        iGrid_left=iGrid;
        bp_remaining = shuffle_bin_radius;
        bp_prev = focal_point;
        total_bp_added = 0;
        while ((0 < bp_remaining) & (0 <= iGrid_left)) {
            bp_to_add = (bp_prev - L_grid(iGrid_left));
            if ((bp_remaining - bp_to_add) < 0) {
                bp_to_add = bp_remaining;
                bp_remaining = 0;
            } else {
                bp_remaining = bp_remaining - bp_to_add;
            }
            // add bit
            smoothed_rate(iGrid) = smoothed_rate(iGrid) + \
                bp_to_add * sigma_rate(iGrid_left);
            // move
            total_bp_added += bp_to_add;            
            bp_prev = L_grid(iGrid_left);
            iGrid_left = iGrid_left - 1;                
        }
        //
        // right
        //
        iGrid_right = iGrid + 1;
        bp_remaining = shuffle_bin_radius;
        bp_prev = focal_point;            
        while ((0 < bp_remaining) & (iGrid_right < nGrids)) {
            bp_to_add = (L_grid(iGrid_right) - bp_prev);
            if ((bp_remaining - bp_to_add) < 0) {
                bp_to_add = bp_remaining;
                bp_remaining = 0;
            } else {
                bp_remaining = bp_remaining - bp_to_add;
            }
            // add bit
            smoothed_rate(iGrid) = smoothed_rate(iGrid) +       \
                bp_to_add * sigma_rate(iGrid_right - 1);
            // move
            total_bp_added += bp_to_add;
            bp_prev = L_grid(iGrid_right);
            iGrid_right = iGrid_right + 1;                
        }
        //
        // normalize
        //
        smoothed_rate(iGrid) /= total_bp_added;
    }
    return(smoothed_rate);
}
