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
Rcpp::NumericVector rcpp_calculate_hwe_p(const Rcpp::IntegerVector reference_hap) {

  Rcpp::IntegerVector x = reference_hap;
  if (x[2] > x[0]) {
    int a = x[2];
    x[2] = x[0];
    x[0] = a;
  }

  int nAA = x[0]; // individuals
  int nAB = x[1]; // individuals
  int nBB = x[2]; // individuals
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

  // determine a mid point
  Rcpp::NumericVector probs = Rcpp::NumericVector(max_het + 1); //note - this is 0-based
  probs[mid] = 1; // made 0-based for C++

  // use recurrence relation - going down  
  int n_het = mid;
  int n_hom_alt = (nBB * 2 + nAB - n_het) / 2;
  int n_hom_ref = (n - n_het - n_hom_alt);
  if ((mid - 2) >= min_het) {
      // am here, fix seq below
      for(int het = mid - 2; het >= min_het; het-=2) {
          probs[het] = probs[het + 2] *                                 \
              double(n_het * (n_het-1)) / double(4 * (n_hom_ref + 1) * (n_hom_alt + 1));
          n_het = n_het - 2;
          n_hom_ref = n_hom_ref + 1;
          n_hom_alt = n_hom_alt + 1;
      }
  }
  
  // use recurrence relationship - going up
  n_het = mid;
  n_hom_alt = (nBB * 2 + nAB - n_het) / 2;
  n_hom_ref = (n - n_het - n_hom_alt);
  if ((mid + 2) <= max_het) {
      for(int het = mid + 2; het <= max_het; het+=2) {            
          probs[het] = probs[het - 2] * \
              (4 * n_hom_ref * n_hom_alt) / ( (n_het + 2) * (n_het + 1));
          n_het = n_het + 2;
          n_hom_ref = n_hom_ref - 1;
          n_hom_alt = n_hom_alt - 1;
      }
  }

  double probs_sum = 0;
  for(int k = 0; k <= max_het; k++) {
      probs_sum = probs_sum + probs[k];
  }
  
  Rcpp::NumericVector all_probs = Rcpp::NumericVector(max_het + 1);
  for(int k = 0; k <= max_het; k++) {
      all_probs[k] = probs[k] / probs_sum;
  }

  double p2 = 0;
  for(int k = 0; k <= max_het; k++) {
      if (all_probs[k] <= all_probs[nAB]) {
          p2 = p2 + all_probs[k];
      }
  }

  // R is being weird about returning straight doubles
  Rcpp::NumericVector p3 = Rcpp::NumericVector(1);
  if (p2 > 1.0) {
      p3 = 1.0;
  } else {
      p3 = p2;
  }
  return(p3);
  
}
