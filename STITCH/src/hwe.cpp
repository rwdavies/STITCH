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
#include <math.h>
using namespace Rcpp;


//' @export
// [[Rcpp::export]]
double rcpp_calculate_hwe_p(const Rcpp::IntegerVector reference_hap) {

  Rcpp::IntegerVector x = reference_hap;
  if (x[2] > x[0]) {
    int a = x[2];
    x[2] = x[0];
    x[0] = a;
  }

  int nAA = x[1]; // individuals
  int nAB = x[2]; // individuals
  int nBB = x[3]; // individuals
  int nA = nAA * 2 + nAB; // number of copies of A allele
  int nB = nBB * 2 + nAB; // number of copies of B allele
  int n = nAA + nAB + nBB; // number of (diploid) individuals

  int min_het = 0;
  int max_het = nAB + 2 * std::min(nAA, nBB);

  int mid = 0;  
  if ( (nA + nB) > 0 ) {
    mid = floor((nA * nB) / ((nA + nB)));
  }

  // make odd if it needs to be
  if ( (nA % 2) != (mid % 2))
    mid = mid + 1;

  

  double p = 1;
  return p;
}
