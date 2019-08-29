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
void Rcpp_rhb_reader_chunk_process(
    arma::imat& rhb,
    arma::imat& hold,
    const Rcpp::StringVector& chunk,
    const int& chunk_length,
    const int& start_snp,
    const int& end_snp,
    Rcpp::IntegerVector& bs,
    Rcpp::IntegerVector& ihold,
    const Rcpp::IntegerVector& haps_to_get,
    const int& final_snp_to_get,
    const int& n_haps,
    const Rcpp::LogicalVector& binary_get_line,
    arma::mat& ref_alleleCount,
    const arma::ivec& rh_in_L,
    Rcpp::LogicalVector& final_snp_gotten
) {
    int k, iiSNP;
    iiSNP = -1;
    while(iiSNP < (chunk_length - 1)) {
        iiSNP++;
        //std::cout << "iiSNP = " << iiSNP << std::endl;
        int iSNP = start_snp - 1 + iiSNP; // 0-based here
        //std::cout << "iSNP = " << iSNP << std::endl;        
        if (binary_get_line(iSNP)) {
            //std::cout << "inside, do hold storage, ihold(0) = " << ihold(0) << std::endl;            
            for(int i_hap = 0; i_hap < n_haps; i_hap++) {
                hold(ihold(0), i_hap) = chunk[iiSNP][2 * haps_to_get(i_hap)] - '0';
            }
            int w = rh_in_L(32 * bs(0) + ihold(0)) - 1;
            ref_alleleCount(w, 0) = sum(hold.row(ihold(0)));
            ref_alleleCount(w, 1) = n_haps;
            ref_alleleCount(w, 2) = ref_alleleCount(w, 0) / ref_alleleCount(w, 1);
            // if full or the final SNP to get, reset bs / ihold
            if ((iSNP == final_snp_to_get) | (ihold(0) == 31)) {
                //std::cout << "do rbs storage" << std::endl;
                // ## overflow or end. note - with reset, should be fine not resetting this
                for(int i_hap = 0; i_hap < n_haps; i_hap++) {
                    std::uint32_t itmp = 0;
                    for (k = 31; k >= 0; k--) {
                        itmp <<= 1;
                        int j = hold(k, i_hap);
                        itmp |= j & 0x1;
                    }
                    rhb(bs(0), i_hap) = itmp; // bs is 0-based
                    // like https://github.com/wch/r-source/blob/5a156a0865362bb8381dcd69ac335f5174a4f60c/src/main/raw.c#L179
                }
                hold.fill(0);
                bs(0)++;
                ihold(0) = 0; // 0-based
                if (iSNP == final_snp_to_get) {
                    iiSNP += 100; // end
                    final_snp_gotten(0) = true;
                }
	    } else {
                ihold(0) += 1;
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
Rcpp::IntegerVector rcpp_int_contract(const arma::ivec& hap) {
    const int nSNPs = hap.size();
    const int nbSNPs = std::ceil(double(nSNPs) / double(32));
    Rcpp::IntegerVector hapc(nbSNPs);
    int imax;
    for(int bs = 0; bs < nbSNPs; bs++) {
        int d32_times_bs = 32 * bs;
        if (bs < (nbSNPs - 1)) {
            imax = 31;
        } else {
            // final one!
            imax = nSNPs - d32_times_bs - 1;
        }
        std::uint32_t itmp = 0;
        //std::cout << "k = ";
        for (int k = imax; k >= 0; k--) {
            // std::cout << k << ", ";
            itmp <<= 1;
            int j = hap(d32_times_bs + k);
            itmp |= j & 0x1;
        }
        //std::cout << std::endl;
        hapc(bs) = itmp;
        // see reader_chunk_process 
    }
    return(hapc);
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


//' @export
// [[Rcpp::export]]
arma::imat inflate_fhb(
    arma::imat& rhb,
    Rcpp::IntegerVector& haps_to_get,
    const int nSNPs
) {
    const int K = rhb.n_cols;
    const int nbSNPs = rhb.n_rows;
    // i think this function might work without kmax
    // but probably safer / simpler to keep it in
    int imax;
    int n_haps_to_get = haps_to_get.size();    
    arma::imat rhi_subset(nSNPs, n_haps_to_get);
    // outer loop is on k, as the columns are "K", i.e. haps
    for(int ik = 0; ik < n_haps_to_get; ik++) {
        int k = haps_to_get(ik);
        for(int bs = 0; bs < nbSNPs; bs++) {
            int d32_times_bs = 32 * bs;
            if (bs < (nbSNPs - 1)) {
                imax = 32;
            } else {
                // final one!
                imax = nSNPs - d32_times_bs;
            }
	    std::uint32_t tmp(rhb(bs, k));
	    // this weird looking code taken largely from R c code
	    // e.g. search for this in R github
	    // SEXP attribute_hidden do_intToBits(SEXP call, SEXP op, SEXP args, SEXP env)
	    for (int i = 0; i < imax; i++, tmp >>= 1) {
	      // might be inefficient
	      // might need to work with temporary matrix?
	      // revisit later if actually slow!
                rhi_subset(d32_times_bs + i, ik) = tmp & 0x1;
	    }
	}
    }
    return(rhi_subset);
}

