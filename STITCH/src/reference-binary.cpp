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
void Rcpp_rhb_reader_chunk_process(
    arma::imat& rhb,
    arma::imat& hold,
    const Rcpp::StringVector& chunk,
    const int& chunk_length,
    const int& start_snp,
    const int& end_snp,
    int& bs,
    int& ihold,
    Rcpp::IntegerVector& haps_to_get,
    int& final_snp_to_get,
    int& n_haps,
    Rcpp::LogicalVector& binary_get_line
) {
    int k, iiSNP;
    iiSNP = -1;
    std::cout << "inside!" << std::endl;
    while(iiSNP < chunk_length) {
        iiSNP++;
        std::cout << "iiSNP = " << iiSNP << std::endl;        
        int iSNP = start_snp - 1 + iiSNP; // 0-based here
        std::cout << "iSNP = " << iSNP << std::endl;
        if (binary_get_line(iSNP)) {
            for(int i_hap = 0; i_hap < n_haps; i_hap++) {
                hold(ihold, i_hap) = chunk[iiSNP][2 * haps_to_get(i_hap)] - '0';
            }
            // if full or the final SNP to get, reset bs / ihold
            if ((iSNP == final_snp_to_get) | (ihold == 31)) {






                
                std::cout << "am here! sort this bad boy ouT!" << std::endl;
                so close!!!!



                
                // ## overflow or end. note - with reset, should be fine not resetting this
                for(int i_hap = 0; i_hap < n_haps; i_hap++) {
                    std::uint32_t itmp = 0;
                    for (k = 31; k >= 0; k--) {
                        itmp <<= 1;
                        int j = hold(k, i_hap);
                        itmp |= j & 0x1;
                    }
                    rhb(bs, k) = itmp; // bs is 0-based
                    // like https://github.com/wch/r-source/blob/5a156a0865362bb8381dcd69ac335f5174a4f60c/src/main/raw.c#L179
                }
                hold.fill(0);
                bs++;
                ihold = 0; // 0-based
                iiSNP += 100; // finish off
	    } else {
                ihold++;
            }
        }
    }
    return;
}



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
