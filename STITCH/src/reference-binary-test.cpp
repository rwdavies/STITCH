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
#include <bitset>


//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_f3_t(arma::mat& rh_numeric_t, arma::rowvec& hap) {
  const int K = rh_numeric_t.n_rows;
  Rcpp::NumericVector x(K);
  for(int k = 0; k < K; k++) {
    x(k) = sum(abs(hap - rh_numeric_t.row(k)));
  }
  return(x);
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_f3(arma::mat& rh_numeric, arma::colvec& hap) {
  const int K = rh_numeric.n_cols;
  Rcpp::NumericVector x(K);
  for(int k = 0; k < K; k++) {
    x(k) = sum(abs(hap - rh_numeric.col(k)));
  }
  return(x);
}





// //' @export
// // [[Rcpp::export]]
// void rcpp_int_expand(arma::ivec& o, arma::ivec& x, const int bSNPs) {
//   int j = 0;
//   for(int bs = 0; bs < bSNPs; bs++) {
//     std::uint32_t tmp(x(bs));
//     for (int k = 0; k < 32; k++, tmp >>= 1) {
//       o(j++) = tmp & 0x1;
//     }
//   }
//   return;
// }

//' @export
// [[Rcpp::export]]
void rcpp_int_expand2(arma::ivec& o, arma::imat& x, const int k, const int bSNPs) {
  int j = 0;
  for(int bs = 0; bs < bSNPs; bs++) {
    std::uint32_t tmp(x(bs, k));
    for (int k = 0; k < 32; k++, tmp >>= 1) {
      o(j++) = tmp & 0x1;
    }
  }
  return;
}


// where rhb is an [nbSNPs x K] matrix

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_f3b(arma::imat& rhb, arma::vec& hap) {
    const int K = rhb.n_cols;
    const int bSNPs = rhb.n_rows;
    Rcpp::NumericVector out(K);
    for(int k = 0; k < K; k++) {
        int j = 0;
	for(int bs = 0; bs < bSNPs; bs++) {
	    std::uint32_t tmp(rhb(bs, k));
	    double d = 0;
	    for (int i = 0; i < 32; i++, tmp >>= 1) {
	        d += std::abs(hap(32 * bs + i) - (tmp & 0x1));
	    }
	    out(k) += d;
	}
    }
    return(out);
}


//' @export
// [[Rcpp::export]]
arma::colvec rcpp_f4b_t(arma::imat& rhb_t, arma::vec& hap) {
    const int K = rhb_t.n_rows;
    const int bSNPs = rhb_t.n_cols;
    arma::colvec out(K);
    out.fill(0);
    // outer loop is on bs / "SNPs"
    for(int bs = 0; bs < bSNPs; bs++) {
        arma::vec h = hap.subvec(32 * bs, 32 * bs + 31);
	arma::ivec rhb_t_colvec = rhb_t.col(bs);
	for(int k = 0; k < K; k++) {
	    std::uint32_t tmp(rhb_t_colvec(k));
	    double d = 0;
	    // this weird looking code taken largely from R c code
	    // e.g. search for this in R github
	    // SEXP attribute_hidden do_intToBits(SEXP call, SEXP op, SEXP args, SEXP env)
	    for (int i = 0; i < 32; i++, tmp >>= 1) {
	        d += std::abs(h(i) - (tmp & 0x1));
	    }
	    out(k) += d;	  
	}
    }
    return(out);
}




//' @export
// [[Rcpp::export]]
arma::rowvec rcpp_f4(arma::mat& rh_numeric, arma::rowvec& hap) {
  const int K = rh_numeric.n_cols;
  const int nSNPs = rh_numeric.n_rows;
  arma::rowvec x(K);
  x.fill(0);
  double h;
  for(int iSNP = 0; iSNP < nSNPs; iSNP++) {
    h = hap(iSNP);
    x += abs(h - rh_numeric.row(iSNP));
  }
  return(x);
}

//' @export
// [[Rcpp::export]]
arma::colvec rcpp_f4_t(arma::mat& rh_numeric_t, arma::colvec& hap) {
  const int nSNPs = rh_numeric_t.n_cols;
  const int K = rh_numeric_t.n_rows;
  arma::colvec x(K);
  x.fill(0);
  double h;
  for(int iSNP = 0; iSNP < nSNPs; iSNP++) {
    h = hap(iSNP);
    x += abs(h - rh_numeric_t.col(iSNP));
  }
  return(x);
}

