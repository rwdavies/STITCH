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


double print_times(
    double prev,
    int suppressOutput,
    std::string past_text,
    std::string next_text
);

arma::mat make_gammaEK_t_from_gammaK_t(
    const arma::mat& gammaK_t,
    const int K,
    const Rcpp::IntegerVector& grid,
    const int snp_start_1_based,
    const int snp_end_1_based,
    double& prev,
    const int suppressOutput,
    std::string& prev_section,
    std::string& next_section,
    const int grid_offset = 0
);



//' @export
// [[Rcpp::export]]
void Rcpp_ref_run_forward_haploid(
    arma::mat& alphaHat_t,
    arma::rowvec& c,
    const arma::mat& eMatGrid_t,
    const arma::cube& alphaMatCurrent_tc,
    const arma::cube& transMatRate_tc_H,
    const arma::mat& priorCurrent_m,
    const int s
) {
    const int K = alphaMatCurrent_tc.n_rows;
    const int nGrids = alphaMatCurrent_tc.n_cols + 1;    
    //
    // initialize
    //
    int k;
    for(k = 0; k < K; k++) {
      alphaHat_t(k, 0) = priorCurrent_m(k, s) * eMatGrid_t(k, 0);
    }
    c(0) = 1 / sum(alphaHat_t.col(0));
    alphaHat_t.col(0) = alphaHat_t.col(0) * c(0);
    //
    for(int iGrid = 1; iGrid < nGrids; iGrid++) {
        // NOTE - previously used this code, which is mathematically right
        // BUT since scaling is being used here, arma::sum(alphaHat_t.col(t - 1) is equal to 1 by definition
        // so can use the below (uncommented) code to simplify
        // alphaConst = transMatRate_t_H(1, t-1) * arma::sum(alphaHat_t.col(t - 1));
        //
        alphaHat_t.col(iGrid) = eMatGrid_t.col(iGrid) % (		   \
            transMatRate_tc_H(0, iGrid - 1, s) * alphaHat_t.col(iGrid - 1) + \
            transMatRate_tc_H(1, iGrid - 1, s) * alphaMatCurrent_tc.slice(s).col(iGrid - 1) );
        c(iGrid) = 1 / arma::sum(alphaHat_t.col(iGrid));
        alphaHat_t.col(iGrid) *= c(iGrid);
    }
    return ;
}

void Rcpp_run_backward_haploid(
    arma::mat& betaHat_t,
    arma::rowvec& c,
    const arma::mat& eMatGrid_t,
    const arma::cube& alphaMatCurrent_tc,
    const arma::cube& transMatRate_tc_H,
    const int s
);




//' @export
// [[Rcpp::export]]
Rcpp::List rcpp_make_sampleReads_from_hap(
    const Rcpp::IntegerVector non_NA_cols,
    const int reference_phred,
    const Rcpp::IntegerVector reference_hap
) {
    Rcpp::List sampleReads(non_NA_cols.length());
    int iSNP;
    for(int i = 0; i < non_NA_cols.length(); i++) {
      iSNP = non_NA_cols[i] - 1; // this is 0-based (note - grid comes later)
        sampleReads[i]=Rcpp::List::create(0, iSNP, reference_phred * (2 * reference_hap[iSNP] - 1), iSNP);
    }
    return sampleReads;
}


// since reference_phred is fixed, alot of these matrices are the same
// pre-calculate them here

//' @export
// [[Rcpp::export]]
void ref_make_ehh(
    const arma::cube& eHapsCurrent_tc,
    const Rcpp::IntegerVector& non_NA_cols,
    arma::cube& ehh_h1_A,
    arma::cube& ehh_h1_S,
    arma::cube& ehh_h0_A,
    arma::cube& ehh_h0_S,
    double reference_phred
) {
    const int S = eHapsCurrent_tc.n_slices;
    int iSNP;
    int nnSNPs = non_NA_cols.length();
    double eps = pow(10,(double(-reference_phred) / 10));
    double pA, pR;
    for(int s = 0; s < S; s++) {
      // do all sites - ignore potential slight over-multiplication, likely faster
      // do h1 first
      pA = 1 - eps;
      pR = eps / 3;
      ehh_h1_A.slice(s) = eHapsCurrent_tc.slice(s) * pA;
      ehh_h1_S.slice(s) = ehh_h1_A.slice(s) + (1 - eHapsCurrent_tc.slice(s)) * pR;
      // do h0
      pA = eps / 3;
      pR = 1 - eps;      
      ehh_h0_A.slice(s) = eHapsCurrent_tc.slice(s) * pA;
      ehh_h0_S.slice(s) = ehh_h0_A.slice(s) + (1 - eHapsCurrent_tc.slice(s)) * pR;      
    }
    return;
}


//' @export
// [[Rcpp::export]]
void rcpp_ref_bound_eMatGrid_t(
    arma::mat& eMatGrid_t,
    const double maxEmissionMatrixDifference,
    bool rescale,
    bool bound
) {
    // do not need to
    if ((1 / maxEmissionMatrixDifference) < (eMatGrid_t.min())) {
        bound = false;
    }
    const int K = eMatGrid_t.n_rows;
    const int nGrids = eMatGrid_t.n_cols;
    double rescale_val, x;
    double d2 = 1 / maxEmissionMatrixDifference;
    double d3 = d2;
    int k;
    for(int iGrid = 0; iGrid < nGrids; iGrid++) {
        if (eMatGrid_t(0, iGrid) < 1) {
            x = eMatGrid_t.col(iGrid).max();
	    rescale_val = 1 / x;	
	    if (rescale) {
	      eMatGrid_t.col(iGrid) *= rescale_val;
	      rescale_val = 1;
	    } else {
	      d3 = rescale_val / d2;	      
	    }
	    if (bound) {
	        for (k = 0; k < K; k++) {
		    if(eMatGrid_t(k, iGrid) < d3) {
		        eMatGrid_t(k, iGrid) = d3;
		    }
		}
	    }
	}
    }
    return;
}

//' @export
// [[Rcpp::export]]
void rcpp_ref_make_eMatGrid_t(
    arma::mat& eMatGrid_t,
    Rcpp::IntegerMatrix& reference_haps,
    const Rcpp::IntegerVector& non_NA_cols,
    const arma::cube& eHapsCurrent_tc,
    const Rcpp::IntegerVector& grid,    
    const int reference_phred,    
    const int s,
    const int iSample,
    const double maxEmissionMatrixDifference,
    arma::cube& ehh_h1_A,
    arma::cube& ehh_h1_S,
    arma::cube& ehh_h0_A,
    arma::cube& ehh_h0_S,
    const bool rescale = true,
    const bool bound = true
) {
    //
    // non_NA_cols is 1-based, which columns are we running over
    //
    const int K = eMatGrid_t.n_rows;
    int iSNP, iGrid, h;
    double eps = pow(10,(double(-reference_phred) / 10));
    double pA, pR, eps2, d2, x, rescale_val;
    int k, cur_grid;
    int prev_grid = -1;
    int nnSNPs = non_NA_cols.length();
    //
    //
    for(int iiSNP = 0; iiSNP < nnSNPs; iiSNP++) {
        iSNP = non_NA_cols(iiSNP) - 1; // here, 0-based
	iGrid = grid(iSNP);
	h = reference_haps(iSNP, iSample); // here, 0 = ref, 1 = alt
	if (h == 1) {
	    eMatGrid_t.col(iGrid) %= ehh_h1_S.slice(s).col(iSNP);
	} else {
            eMatGrid_t.col(iGrid) %= ehh_h0_S.slice(s).col(iSNP);
	}
	//std::cout << "iSNP = " << iSNP << ", iGrid = " << iGrid << ", h = " << h << std::endl;
	// if (h == 1) {
	//   pA = 1 - eps;
	//   pR = eps / 3;
	// } else {
	//   pA = eps / 3;
	//   pR = 1 - eps;
	// }
	//eMatGrid_t.col(iGrid) %= (eHapsCurrent_tc.slice(s).col(iSNP) * pA + (1 - eHapsCurrent_tc.slice(s).col(iSNP)) * pR);
	// eHapsHelper_A.col(iSNP) = eHapsCurrent_tc.slice(s).col(iSNP) * pA;
	// eHapsHelper_B.col(iSNP) = (1 - eHapsCurrent_tc.slice(s).col(iSNP)) * pR;
	// //
	// eMatGrid_t.col(iGrid) %= (eHapsHelper_A.col(iSNP) + eHapsHelper_B.col(iSNP));
	// 
	// otherwise, not too worried
	// cur_grid = iGrid;
	// if (cur_grid == prev_grid) {
	//   eps2 *= eps;
	// } else {
	//   eps2 = 1;
	// }
	// prev_grid = cur_grid;
	// screw it, be inefficient!
    }
    if (rescale | bound) {
      rcpp_ref_bound_eMatGrid_t(eMatGrid_t, maxEmissionMatrixDifference, rescale, bound);
    }
    return;
}



//' @export
// [[Rcpp::export]]
void ref_make_haploid_gammaUpdate_t(
    int s,
    arma::cube& gammaSum0_tc,
    arma::cube& gammaSum1_tc,    
    const arma::mat& gamma_t,
    const arma::cube& eHapsCurrent_tc,
    const Rcpp::IntegerMatrix& reference_haps,
    const Rcpp::IntegerVector& non_NA_cols,
    const int iSample,
    const Rcpp::IntegerVector& grid,
    const double reference_phred,
    arma::cube& ehh_h1_A,
    arma::cube& ehh_h1_S,
    arma::cube& ehh_h0_A,
    arma::cube& ehh_h0_S
) {
    //
    //const int nSNPs = eHapsCurrent_t.n_cols;
    //
    arma::ivec bqU, pRU;
    arma::colvec gamma_t_col;
    //double d3, eps, a1, a2, y, d, d1, d2, b;
    double eps = pow(10,(double(-reference_phred) / 10));
    double pR = -1;
    double pA = -1;
    int iGrid_prev = -1;
    int h, iGrid, iSNP;
    arma::colvec a1_col, a2_col, y_col, d_col;
    arma::colvec b_col, d1_col, d2_col;
    int nnSNPs = non_NA_cols.length();    
    //
    // here, only need to loop over non-empty eMatGrid_t
    //
    for(int iiSNP = 0; iiSNP < nnSNPs; iiSNP++) {
        iSNP = non_NA_cols(iiSNP) - 1; // here, 0-based
	iGrid = grid(iSNP);
        if (iGrid > iGrid_prev) {
            // do not always need to update
            gamma_t_col = gamma_t.col(iGrid);
        }
        iGrid_prev = iGrid;	
	h = reference_haps(iSNP, iSample); // here, 0 = ref, 1 = alt
	// if (h == 1) {
	//   pA = 1 - eps;
	//   pR = eps / 3;
	// } else {
	//   pA = eps / 3;
	//   pR = 1 - eps;
	// }
	//a1_col = pA * eHapsCurrent_tc.slice(s).col(iSNP);
	//a2_col = pR * (1 - eHapsCurrent_tc.slice(s).col(iSNP));
	if (h == 1) {
	  a1_col = ehh_h1_A.slice(s).col(iSNP);
	  y_col  = ehh_h1_S.slice(s).col(iSNP);
	} else {
	  a1_col = ehh_h0_A.slice(s).col(iSNP);
	  y_col  = ehh_h0_S.slice(s).col(iSNP);
	}
	// a1_col = eHapsHelper_A.col(iSNP);
	// a2_col = eHapsHelper_B.col(iSNP);
	// y_col = a1_col + a2_col;
	d_col = gamma_t_col / y_col;
	gammaSum0_tc.slice(s).col(iSNP) += a1_col % d_col; // despite being "gamma", this is per-SNP
	gammaSum1_tc.slice(s).col(iSNP) += y_col % d_col;
    }
    return;
}





// forwardBackwardHaploid for reference
// different input
// streamlined
// only needs to update
// and maybe output gamma (?efficiency)

//' @export
// [[Rcpp::export]]
Rcpp::List reference_fbh(
    const arma::cube& eHapsCurrent_tc,
    const arma::cube& alphaMatCurrent_tc,
    const arma::cube& transMatRate_tc_H,
    const arma::mat& priorCurrent_m,
    arma::mat& alphaHat_t,
    arma::mat& betaHat_t,
    arma::mat& gamma_t,
    arma::mat& eMatGrid_t,
    arma::cube& ehh_h1_A,
    arma::cube& ehh_h1_S,
    arma::cube& ehh_h0_A,
    arma::cube& ehh_h0_S,
    const double maxDifferenceBetweenReads,
    const double maxEmissionMatrixDifference,
    const int suppressOutput,
    Rcpp::IntegerMatrix& reference_haps,
    const Rcpp::IntegerVector& non_NA_cols,
    const int iSample,
    const double reference_phred,
    arma::cube& gammaSum0_tc,
    arma::cube& gammaSum1_tc,    
    arma::cube& alphaMatSum_tc,
    arma::cube& hapSum_tc,
    arma::mat& priorSum_m,
    const bool return_extra = false,
    const bool return_gammaK = false,    
    const Rcpp::IntegerVector grid = 0,
    const bool rescale_eMatGrid_t = true,
    const bool bound_eMatGrid_t = true,
    const bool run_fb_subset = false    
) {
  double prev=clock();
  std::string prev_section="Null";
  std::string next_section="Initialize variables";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  //
  // constants
  //
  //
  const int nGrids = alphaMatCurrent_tc.n_cols + 1;  // what we iterate over / grid
  const int K = eHapsCurrent_tc.n_rows; // traditional K for haplotypes
  const int S = eHapsCurrent_tc.n_slices;
  const int nSNPs = eHapsCurrent_tc.n_cols;
  //
  // new variables
  //
  arma::rowvec c = arma::zeros(1, nGrids);
  // variables for transition matrix and initialization
  // int variables and such
  int iGrid, k, s;
  Rcpp::List to_return;
  Rcpp::List list_of_gammaK_t;  
  //  
  Rcpp::NumericVector alphaStart, betaEnd;
  arma::mat alphaHatBlocks_t, betaHatBlocks_t;
  double g_temp = 0;
  //
  // everything works on s here
  //
  for(s = 0; s < S; s++) {
      // always reset - always passed in
      alphaHat_t.fill(0);
      betaHat_t.fill(0);
      eMatGrid_t.fill(1);
      //
      // so only need eMatGrid_t marginally
      //
      next_section="Make eMatGrid_t for special reference";
      prev=print_times(prev, suppressOutput, prev_section, next_section);
      prev_section=next_section;
      rcpp_ref_make_eMatGrid_t(eMatGrid_t, reference_haps, non_NA_cols, eHapsCurrent_tc, grid, reference_phred, s, iSample, maxEmissionMatrixDifference, ehh_h1_A, ehh_h1_S, ehh_h0_A, ehh_h0_S, rescale_eMatGrid_t, bound_eMatGrid_t);
      //
      //
      next_section="Forward recursion";
      prev=print_times(prev, suppressOutput, prev_section, next_section);
      prev_section=next_section;
      //
      Rcpp_ref_run_forward_haploid(alphaHat_t, c, eMatGrid_t, alphaMatCurrent_tc, transMatRate_tc_H, priorCurrent_m, s);
      //
      // backward recursion
      //
      next_section="Backward recursion";
      prev=print_times(prev, suppressOutput, prev_section, next_section);
      prev_section=next_section;
      betaHat_t.col(nGrids - 1).fill(c(nGrids - 1));
      Rcpp_run_backward_haploid(betaHat_t, c, eMatGrid_t, alphaMatCurrent_tc, transMatRate_tc_H, s);
      //
      // make gamma
      //
      next_section="Make gamma";
      prev=print_times(prev, suppressOutput, prev_section, next_section);
      prev_section=next_section;
      //
      gamma_t = alphaHat_t % betaHat_t;
      // normalize as well
      for(iGrid = 0; iGrid < nGrids; iGrid++) {
          g_temp = 1 / c(iGrid);
          gamma_t.col(iGrid) *= g_temp;
      }
      //
      // inline - make changes
      //
      next_section="Update priorSum";
      prev=print_times(prev, suppressOutput, prev_section, next_section);
      prev_section=next_section;
      //
      priorSum_m.col(s) += gamma_t.col(0);
      //
      // update hapSum
      //
      next_section="Update hapSum";
      prev=print_times(prev, suppressOutput, prev_section, next_section);
      prev_section=next_section;
      //
      hapSum_tc.slice(s) += gamma_t;
      //
      // update gammaUpdate
      //
      next_section="Update gammaUpdate";
      prev=print_times(prev, suppressOutput, prev_section, next_section);
      prev_section=next_section;
      //
      ref_make_haploid_gammaUpdate_t(s, gammaSum0_tc, gammaSum1_tc, gamma_t, eHapsCurrent_tc, reference_haps, non_NA_cols, iSample, grid, reference_phred, ehh_h1_A, ehh_h1_S, ehh_h0_A, ehh_h0_S);
      //
      // update alphaMatSum
      //
      next_section="Update alphaMatSum";
      prev=print_times(prev, suppressOutput, prev_section, next_section);
      prev_section=next_section;
      for(int iGrid = 0; iGrid < nGrids - 1; iGrid++) {
          alphaMatSum_tc.slice(s).col(iGrid) += transMatRate_tc_H(1, iGrid, s) * (alphaMatCurrent_tc.slice(s).col(iGrid) % betaHat_t.col(iGrid + 1) % eMatGrid_t.col(iGrid + 1));
      }
      //
      if (return_gammaK) {
          list_of_gammaK_t.push_back(gamma_t);
      }
  }
  next_section="Finalize";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  if (return_gammaK) {
    to_return.push_back(list_of_gammaK_t, "list_of_gammaK_t");
  }
  if (return_extra) {
    to_return.push_back(eMatGrid_t, "eMatGrid_t");
  }
  //
  //
  next_section="Done";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  return(to_return);  
}
