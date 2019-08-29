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
    const arma::cube& alphaMatCurrentX_tc,
    const arma::cube& transMatRate_tc_H,
    const arma::mat& priorCurrent_m,
    const int s
) {
    const int K = alphaMatCurrentX_tc.n_rows;
    const int nGrids = alphaMatCurrentX_tc.n_cols + 1;    
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
            alphaMatCurrentX_tc.slice(s).col(iGrid - 1) );
	//            transMatRate_tc_H(1, iGrid - 1, s) * alphaMatCurrent_tc.slice(s).col(iGrid - 1) );	
        c(iGrid) = 1 / arma::sum(alphaHat_t.col(iGrid));
        alphaHat_t.col(iGrid) *= c(iGrid);
    }
    return ;
}



//' @export
// [[Rcpp::export]]
void Rcpp_ref_run_backward_haploid(
    arma::mat& betaHat_t,
    arma::rowvec& c,
    const arma::mat& eMatGrid_t,
    const arma::cube& alphaMatCurrentX_tc,
    const arma::cube& transMatRate_tc_H,
    arma::cube& alphaMatSum_tc,
    const int s
) {
    const int nGrids = eMatGrid_t.n_cols;
    double x;
    arma::colvec e_times_b;
    for(int iGrid = nGrids - 2; iGrid >= 0; --iGrid) {
        // now using alphaMatCurrentX which has transMat included with it
        //x = transMatRate_tc_H(1, iGrid, s) * sum(alphaMatCurrent_tc.slice(s).col(iGrid) % e_times_b);	
        e_times_b = eMatGrid_t.col(iGrid + 1) % betaHat_t.col(iGrid + 1);
        x = sum(alphaMatCurrentX_tc.slice(s).col(iGrid) % e_times_b);
        betaHat_t.col(iGrid) = c(iGrid) * (x + transMatRate_tc_H(0, iGrid, s) * e_times_b);
	//
	alphaMatSum_tc.slice(s).col(iGrid) += e_times_b;
	// now - alphaMatSum update included here!
    }
    return;
}




//' @export
// [[Rcpp::export]]
Rcpp::List rcpp_make_sampleReads_from_hap(
    const Rcpp::IntegerVector rh_in_L,
    const int reference_phred,
    const Rcpp::IntegerVector reference_hap
) {
    Rcpp::List sampleReads(rh_in_L.length());
    int iSNP;
    for(int i = 0; i < rh_in_L.length(); i++) {
        iSNP = rh_in_L[i] - 1; // this is 0-based (note - grid comes later)
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
    arma::cube& ehh_h1_S,
    arma::cube& ehh_h1_D,
    arma::cube& ehh_h0_S,
    arma::cube& ehh_h0_D,
    double reference_phred
) {
    const int S = eHapsCurrent_tc.n_slices;
    double eps = pow(10,(double(-reference_phred) / 10));
    double pA, pR;
    for(int s = 0; s < S; s++) {
      // do all sites - ignore potential slight over-multiplication, likely faster
      // do h1 first
      pA = 1 - eps;
      pR = eps / 3;
      ehh_h1_S.slice(s) = eHapsCurrent_tc.slice(s) * pA + (1 - eHapsCurrent_tc.slice(s)) * pR;
      ehh_h1_D.slice(s) = (eHapsCurrent_tc.slice(s) * pA) / ehh_h1_S.slice(s);
      // do h0
      pA = eps / 3;
      pR = 1 - eps;      
      ehh_h0_S.slice(s) = eHapsCurrent_tc.slice(s) * pA + (1 - eHapsCurrent_tc.slice(s)) * pR;
      ehh_h0_D.slice(s) = (eHapsCurrent_tc.slice(s) * pA) / ehh_h0_S.slice(s);
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
    Rcpp::IntegerVector& reference_hap,
    const Rcpp::IntegerVector& rh_in_L,
    const arma::cube& eHapsCurrent_tc,
    const Rcpp::IntegerVector& grid,    
    const int reference_phred,    
    const int s,
    const double maxEmissionMatrixDifference,
    arma::cube& ehh_h1_S,
    arma::cube& ehh_h0_S,
    const bool rescale = true,
    const bool bound = true
) {
    //
    // rh_in_L is 1-based, which columns are we running over
    //
    const int K = eMatGrid_t.n_rows;
    int iSNP, iGrid, h;
    // double eps = pow(10,(double(-reference_phred) / 10));
    //int k, cur_grid;
    //int prev_grid = -1;
    int nnSNPs = rh_in_L.length();
    //
    //
    double d2 = 1 / maxEmissionMatrixDifference;
    double d3 = d2;
    double rescale_val;
    //    
    for(int iiSNP = 0; iiSNP < nnSNPs; iiSNP++) {
        iSNP = rh_in_L(iiSNP) - 1; // here, 0-based
	iGrid = grid(iSNP);
	h = reference_hap(iiSNP); // here, 0 = ref, 1 = alt
	if (h == 1) {
	    eMatGrid_t.col(iGrid) %= ehh_h1_S.slice(s).col(iSNP);
	} else {
            eMatGrid_t.col(iGrid) %= ehh_h0_S.slice(s).col(iSNP);
	}
        //
        // now - if this has grid, can overflow. check periodically
        //
        if ((rescale | bound) & (iiSNP > 10)) {
            if (iGrid == grid(rh_in_L(iiSNP - 9) - 1)) {
                // copy
                double x = eMatGrid_t.col(iGrid).max();
                rescale_val = 1 / x;	
                if (rescale) {
                    eMatGrid_t.col(iGrid) *= rescale_val;
                    rescale_val = 1;
                } else {
                    d3 = rescale_val / d2;	      
                }
                if (bound) {
                    for (int k = 0; k < K; k++) {
                        if(eMatGrid_t(k, iGrid) < d3) {
                            eMatGrid_t(k, iGrid) = d3;
                        }
                    }
                }                
            }
        }
        //
        // end of period checker
        //
    }
    if (rescale | bound) {
      rcpp_ref_bound_eMatGrid_t(eMatGrid_t, maxEmissionMatrixDifference, rescale, bound);
    }
    return;
    //
    // old code
    // 
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



//' @export
// [[Rcpp::export]]
void ref_make_haploid_gammaUpdate_t(
    int s,
    arma::cube& gammaSum0_tc,
    arma::cube& gammaSum1_tc,    
    const arma::mat& gamma_t,
    const arma::cube& eHapsCurrent_tc,
    const Rcpp::IntegerVector& reference_hap,
    const Rcpp::IntegerVector& rh_in_L,
    const Rcpp::IntegerVector& grid,
    const double reference_phred,
    arma::cube& ehh_h1_D,    
    arma::cube& ehh_h0_D
) {
    //
    //const int nSNPs = eHapsCurrent_t.n_cols;
    //
    arma::ivec bqU, pRU;
    //double d3, eps, a1, a2, y, d, d1, d2, b;
    // double eps = pow(10,(double(-reference_phred) / 10));
    //double pR = -1;
    //double pA = -1;
    int iGrid_prev = -1;
    int h, iGrid, iSNP;
    //arma::colvec a1_col, a2_col, y_col, d_col, b_col, d1_col, d2_col, gamma_t_col;
    arma::colvec gamma_t_col, d_col;
    int nnSNPs = rh_in_L.length();    
    //
    // here, only need to loop over non-empty eMatGrid_t
    //
    for(int iiSNP = 0; iiSNP < nnSNPs; iiSNP++) {
        iSNP = rh_in_L(iiSNP) - 1; // here, 0-based
	iGrid = grid(iSNP);
        if (iGrid > iGrid_prev) {
            // do not always need to update
            gamma_t_col = gamma_t.col(iGrid);
        }
        iGrid_prev = iGrid;	
	h = reference_hap(iiSNP); // here, 0 = ref, 1 = alt
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
	  d_col  = ehh_h1_D.slice(s).col(iSNP);	  
	} else {
	  d_col  = ehh_h0_D.slice(s).col(iSNP);	  	  
	}
	// a1_col = eHapsHelper_A.col(iSNP);
	// a2_col = eHapsHelper_B.col(iSNP);
	// y_col = a1_col + a2_col;
	gammaSum0_tc.slice(s).col(iSNP) += gamma_t_col % d_col; // despite being "gamma", this is per-SNP
	gammaSum1_tc.slice(s).col(iSNP) += gamma_t_col;
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
    const arma::cube& alphaMatCurrentX_tc,
    const arma::cube& transMatRate_tc_H,
    const arma::mat& priorCurrent_m,
    arma::mat& alphaHat_t,
    arma::mat& betaHat_t,
    arma::mat& gamma_t,
    arma::mat& eMatGrid_t,
    arma::cube& ehh_h1_S,
    arma::cube& ehh_h1_D,
    arma::cube& ehh_h0_S,
    arma::cube& ehh_h0_D,
    const double maxDifferenceBetweenReads,
    const double maxEmissionMatrixDifference,
    int suppressOutput,
    Rcpp::IntegerVector& reference_hap,
    const Rcpp::IntegerVector& rh_in_L,
    const double reference_phred,
    arma::cube& gammaSum0_tc,
    arma::cube& gammaSum1_tc,    
    arma::cube& alphaMatSum_tc,
    arma::cube& hapSum_tc,
    arma::mat& priorSum_m,
    const Rcpp::List& list_of_break_results,
    Rcpp::List& list_of_fromMat,
    Rcpp::List& list_of_fbd_store,
    const Rcpp::IntegerVector& nbreaks,
    const bool save_fbd_store,
    const int iSample1,
    const int pshaM,
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
  const int nGrids = alphaMatCurrentX_tc.n_cols + 1;  // what we iterate over / grid
  //const int K = eHapsCurrent_tc.n_rows; // traditional K for haplotypes
  const int S = eHapsCurrent_tc.n_slices;
  //const int nSNPs = eHapsCurrent_tc.n_cols;
  //
  // new variables
  //
  arma::rowvec c = arma::zeros(1, nGrids);
  // variables for transition matrix and initialization
  // int variables and such
  int iGrid, s;
  Rcpp::List to_return;
  Rcpp::List list_of_gammaK_t;  
  //  
  double g_temp = 0;
  //
  // do inside this
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
      rcpp_ref_make_eMatGrid_t(eMatGrid_t, reference_hap, rh_in_L, eHapsCurrent_tc, grid, reference_phred, s, maxEmissionMatrixDifference, ehh_h1_S, ehh_h0_S, rescale_eMatGrid_t, bound_eMatGrid_t);
      //
      //
      next_section="Forward recursion";
      prev=print_times(prev, suppressOutput, prev_section, next_section);
      prev_section=next_section;
      //
      Rcpp_ref_run_forward_haploid(alphaHat_t, c, eMatGrid_t, alphaMatCurrentX_tc, transMatRate_tc_H, priorCurrent_m, s);
      //
      // backward recursion
      //
      next_section="Backward recursion";
      prev=print_times(prev, suppressOutput, prev_section, next_section);
      prev_section=next_section;
      betaHat_t.col(nGrids - 1).fill(c(nGrids - 1));
      Rcpp_ref_run_backward_haploid(betaHat_t, c, eMatGrid_t, alphaMatCurrentX_tc, transMatRate_tc_H, alphaMatSum_tc, s);
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
      ref_make_haploid_gammaUpdate_t(s, gammaSum0_tc, gammaSum1_tc, gamma_t, eHapsCurrent_tc, reference_hap, rh_in_L, grid, reference_phred, ehh_h1_D, ehh_h0_D);
      //
      // update alphaMatSum
      //
      // NOTE - this is GONE from original version (below)
      // the constant transMatRate multiplication and alphaMatCurrent are outside this loop
      // also, the betaHat times eMatGrid is done when those are multiplied in the backward loop
      //
      // for(int iGrid = 0; iGrid < nGrids - 1; iGrid++) {
      //     alphaMatSum_tc.slice(s).col(iGrid) += \
      // 	    transMatRate_tc_H(1, iGrid, s) *	\
      // 	    (
      // 	     alphaMatCurrent_tc.slice(s).col(iGrid) %			\
      // 	     betaHat_t.col(iGrid + 1) % eMatGrid_t.col(iGrid + 1)	\
      // 	     );
      // }
      //
      if (return_gammaK) {
	next_section="Add to list of gammak_t";
	prev=print_times(prev, suppressOutput, prev_section, next_section);
	prev_section=next_section;
          list_of_gammaK_t.push_back(gamma_t);
      }
      // previous section in R only
      if (0 < nbreaks(s)) {
	next_section="Do fromMat stuff";
	prev=print_times(prev, suppressOutput, prev_section, next_section);
	prev_section=next_section;
	  Rcpp::IntegerMatrix break_results = as<Rcpp::IntegerMatrix>(list_of_break_results[s]);
	  for(int iBreak = 0; iBreak < (nbreaks(s)); iBreak++) {
	    int from = break_results(iBreak, 0); // "left_grid_break_0_based"
	    int to = break_results(iBreak, 3);//"right_grid_break_0_based"
	    arma::colvec hp1 = gamma_t.col(from);
	    arma::rowvec hp2 = gamma_t.col(to).t();
	    arma::cube list_of_fromMat_cube = as<arma::cube>(list_of_fromMat(s));
	    list_of_fromMat_cube.row(iBreak) += hp1 * hp2;
	  }
	  // also, maybe, save fbd_store
	  if (save_fbd_store & (iSample1 <= pshaM)) {
	    Rcpp::List list_of_fbd_store_s = as<Rcpp::List>(list_of_fbd_store(s)); // ugh
	    Rcpp::List temp_list;
	    temp_list.push_back(gamma_t, "gammaK_t");
	    list_of_fbd_store_s(iSample1 - 1) = temp_list;
	  }
      }
  }
  next_section="Finalize and exit";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  if (return_gammaK) {
    to_return.push_back(list_of_gammaK_t, "list_of_gammaK_t");
  }
  if (return_extra) {
    to_return.push_back(eMatGrid_t, "eMatGrid_t");
  }
  return(to_return);  
}


//' @export
// [[Rcpp::export]]
void rcpp_make_alphaMatSumX_tc(
    arma::cube& alphaMatCurrent_tc,
    arma::cube& alphaMatCurrentX_tc,    
    arma::cube& transMatRate_tc_H
) {
    const int S = alphaMatCurrent_tc.n_slices;
    const int nGrids = alphaMatCurrent_tc.n_cols + 1; 
    for(int s = 0; s < S; s++) {
        alphaMatCurrentX_tc.slice(s) = alphaMatCurrent_tc.slice(s);
	for(int iGrid = 0; iGrid < (nGrids - 1); iGrid++) {
	   alphaMatCurrentX_tc.slice(s).col(iGrid) *= transMatRate_tc_H(1, iGrid, s);
	}
    }
    return;
}


//' @export
// [[Rcpp::export]]
void rcpp_finalize_alphaMatSum_tc(
    arma::cube& alphaMatSum_tc,
    arma::cube& transMatRate_tc_H,
    arma::cube& alphaMatCurrentX_tc
) {
    const int S = alphaMatSum_tc.n_slices;
    //const int nGrids = alphaMatSum_tc.n_cols + 1; 
    for(int s = 0; s < S; s++) {
        alphaMatSum_tc.slice(s) %= alphaMatCurrentX_tc.slice(s);
	// for(int iGrid = 0; iGrid < (nGrids - 1); iGrid++) {
	//     alphaMatSum_tc.slice(s).col(iGrid) *= transMatRate_tc_H(1, iGrid, s);
	// }
    }
    return;
}


