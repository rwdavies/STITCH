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
    const Rcpp::NumericVector & sigmaSum,
    const Rcpp::NumericVector & sigma_rate,
    const Rcpp::IntegerVector & L_grid,
    const Rcpp::IntegerVector & grid_distances,
    const int nGrids,
    const int shuffle_bin_radius
) {

    int iSNP;
    int focal_point;
    Rcpp::NumericVector smoothed_rate;
    for(iSNP = 0; iSNP < (nGrids - 1); iSNP++) {
      focal_point = ((L_grid(iSNP) + L_grid(iSNP + 1)) / 2); // automatically rounded
      if (
	  (Rcpp::min(L_grid) <= (focal_point - shuffle_bin_radius)) &
	  ((focal_point + shuffle_bin_radius) <= Rcpp::max(L_grid))
	  ) {
	int a = 1;
      }
    }
    
    
    return(smoothed_rate);
  
}
