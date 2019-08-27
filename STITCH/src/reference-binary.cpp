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


// does not necessarily need to be the most efficient
// start with getting logic right
// later, when it matters, engineer efficient solutions

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector rcpp_int_expand(arma::ivec& hapc, const int nSNPs) {
  const int nbSNPs = hapc.size();
  //const int nSNPs = nbSNPs * 32;
  Rcpp::IntegerVector hap(nSNPs);
  int j = 0;
  int imax;
  for(int bs = 0; bs < nbSNPs; bs++) {
    if (bs < (nbSNPs - 1)) {
      imax = 32;
    } else {
      // final one!
      imax = nSNPs - 32 * bs;
    }
    std::uint32_t tmp(hapc(bs));
    for (int i = 0; i < imax; i++, tmp >>= 1) {
      hap(j++) = tmp & 0x1;
    }
  }
  return(hap);
}


//' @export
// [[Rcpp::export]]
arma::colvec calc_dist_between_rhb_t_and_hap(
    arma::imat& rhb_t,
    arma::vec& hap,
    const int nSNPs
) {
    const int K = rhb_t.n_rows;
    const int nbSNPs = rhb_t.n_cols;
    arma::colvec out(K);
    out.fill(0);
    // i think this function might work without kmax
    // but probably safer / simpler to keep it in
    int imax; 
    // outer loop is on bs / "SNPs"
    for(int bs = 0; bs < nbSNPs; bs++) {
	if (bs < (nbSNPs - 1)) {
	  imax = 32;
	} else {
	  // final one!
	  imax = nSNPs - 32 * bs;
	}
        arma::vec h = hap.subvec(32 * bs, 32 * bs + imax - 1);
	arma::ivec rhb_t_colvec = rhb_t.col(bs);
	for(int k = 0; k < K; k++) {
	    std::uint32_t tmp(rhb_t_colvec(k));
	    double d = 0;
	    // this weird looking code taken largely from R c code
	    // e.g. search for this in R github
	    // SEXP attribute_hidden do_intToBits(SEXP call, SEXP op, SEXP args, SEXP env)
	    for (int i = 0; i < imax; i++, tmp >>= 1) {
	        d += std::abs(h(i) - (tmp & 0x1));
	    }
	    out(k) += d;	  
	}
    }
    return(out);
}



//' @export
// [[Rcpp::export]]
arma::imat inflate_fhb_t(
    arma::imat& rhb_t,
    Rcpp::IntegerVector& haps_to_get,
    const int nSNPs
) {
    const int K = rhb_t.n_rows;
    const int nbSNPs = rhb_t.n_cols;
    // i think this function might work without kmax
    // but probably safer / simpler to keep it in
    int imax;
    int n_haps_to_get = haps_to_get.size();    
    arma::imat rhi_t_subset(n_haps_to_get, nSNPs);
    // outer loop is on bs / "SNPs"
    for(int bs = 0; bs < nbSNPs; bs++) {
	if (bs < (nbSNPs - 1)) {
	  imax = 32;
	} else {
	  // final one!
	  imax = nSNPs - 32 * bs;
	}
	for(int ik = 0; ik < n_haps_to_get; ik++) {
  	    int k = haps_to_get(ik);
	    std::uint32_t tmp(rhb_t(k, bs));
	    // this weird looking code taken largely from R c code
	    // e.g. search for this in R github
	    // SEXP attribute_hidden do_intToBits(SEXP call, SEXP op, SEXP args, SEXP env)
	    for (int i = 0; i < imax; i++, tmp >>= 1) {
	      // might be inefficient
	      // might need to work with temporary matrix?
	      // revisit later if actually slow!
              rhi_t_subset(ik, 32 * bs + i) = tmp & 0x1;
	    }
	}
    }
    return(rhi_t_subset);
}
