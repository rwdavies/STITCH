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

void rcpp_make_eMatRead_t(
    arma::mat& eMatRead_t,
    const Rcpp::List& sampleReads,
    const arma::cube& eHapsCurrent_tc,
    const int s,
    const double maxDifferenceBetweenReads,
    const int Jmax,
    arma::mat& eMatHapOri_t,
    const arma::vec& pRgivenH1,
    const arma::vec& pRgivenH2,
    double& prev,
    int suppressOutput,
    std::string& prev_section,
    std::string& next_section,
    const bool run_pseudo_haploid = false,
    const bool rescale_eMatRead_t = true
);



void Rcpp_run_forward_haploid(
    arma::mat& alphaHat_t,
    arma::rowvec& c,
    const arma::mat& eMatGrid_t,
    const arma::cube& alphaMatCurrent_tc,
    const arma::cube& transMatRate_tc_H,
    const arma::mat& priorCurrent_m,
    const int s,
    const Rcpp::NumericVector alphaStart = 0,
    bool run_fb_subset = false
);

void Rcpp_run_backward_haploid(
    arma::mat& betaHat_t,
    arma::rowvec& c,
    const arma::mat& eMatGrid_t,
    const arma::cube& alphaMatCurrent_tc,
    const arma::cube& transMatRate_tc_H,
    const int s
);

void make_haploid_gammaUpdate_t(
    int s,
    arma::cube& gammaSum0_tc,
    arma::cube& gammaSum1_tc,    
    const Rcpp::List& sampleReads,
    const arma::mat& gamma_t,
    const arma::cube& eHapsCurrent_tc,
    const arma::mat& eMatRead_t,    
    const arma::mat& eMatHapOri_t,
    const arma::vec& pRgivenH1,
    const arma::vec& pRgivenH2,
    const bool run_pseudo_haploid = false
);


void rcpp_make_eMatGrid_t(
    arma::mat& eMatGrid_t,
    const arma::mat& eMatRead_t,
    const Rcpp::IntegerVector& H,
    const Rcpp::List sampleReads,
    const int hap,
    const int nGrids,
    double& prev,
    int suppressOutput,
    std::string& prev_section,
    std::string& next_section,
    const int run_fb_grid_offset = 0,
    const bool use_all_reads = false,
    const bool bound = false,
    const double maxEmissionMatrixDifference = 1000,
    const bool rescale = false
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
    //
    for(int iiSNP = 0; iiSNP < nnSNPs; iiSNP++) {
        iSNP = non_NA_cols(iiSNP) - 1; // here, 0-based
	iGrid = grid(iSNP);
	h = reference_haps(iSNP, iSample); // here, 0 = ref, 1 = alt
	//std::cout << "iSNP = " << iSNP << ", iGrid = " << iGrid << ", h = " << h << std::endl;
	if (h == 1) {
	  pA = 1 - eps;
	  pR = eps / 3;
	} else {
	  pA = eps / 3;
	  pR = 1 - eps;
	}
	eMatGrid_t.col(iGrid) %= ( eHapsCurrent_tc.slice(s).col(iSNP) * pA + (1 - eHapsCurrent_tc.slice(s).col(iSNP)) * pR);
	// 
	// otherwise, not too worried
	cur_grid = iGrid;
	if (cur_grid == prev_grid) {
	  eps2 *= eps;
	} else {
	  eps2 = 1;
	}
	prev_grid = cur_grid;
	// screw it, be inefficient!
	if (
	    (rescale | bound)
	) {
	   eps2 = 1;
	   x = eMatGrid_t.col(iGrid).max();
	   rescale_val = 1 / x;
	   if (rescale) {
	     eMatGrid_t.col(iGrid) *= rescale_val;
	   }
	   if (bound) {
	     d2 = 1 / maxEmissionMatrixDifference;
	     for (k = 0; k < K; k++) {
	       if(eMatGrid_t(k, iGrid) < (d2)) {
		 eMatGrid_t(k, iGrid) = d2;
	       }
	     }
	   }
	}
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
    const Rcpp::List& sampleReads,
    const arma::cube& eHapsCurrent_tc,
    const arma::cube& alphaMatCurrent_tc,
    const arma::cube& transMatRate_tc_H,
    const arma::mat& priorCurrent_m,
    arma::mat& alphaHat_t,
    arma::mat& betaHat_t,
    arma::mat& gamma_t,
    arma::mat& eMatGrid_t,
    const double maxDifferenceBetweenReads,
    const double maxEmissionMatrixDifference,
    const int Jmax,
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
    const bool run_pseudo_haploid,
    const arma::mat& blocks_for_output,
    const Rcpp::List& prev_list_of_alphaBetaBlocks,
    const int i_snp_block_for_alpha_beta = 0,
    const bool generate_fb_snp_offsets = false,
    const bool run_fb_subset = false,
    const int run_fb_grid_offset = 0, // this is 0-based
    const bool return_extra = false,
    const bool return_gamma = false,
    const bool return_gammaK = false,    
    const bool return_hapDosage = true,
    const bool update_in_place = false,
    const bool pass_in_alphaBeta = false, // whether to pass in pre-made alphaHat, betaHat
    const bool output_haplotype_dosages = false, // whether to output state probabilities
    int snp_start_1_based = -1,
    int snp_end_1_based = -1,
    const Rcpp::IntegerVector grid = 0,
    const bool rescale_eMatGrid_t = true,
    const bool bound_eMatGrid_t = true
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
  const int nReads = sampleReads.size();
  const int nGrids = alphaMatCurrent_tc.n_cols + 1;  // what we iterate over / grid
  const int K = eHapsCurrent_tc.n_rows; // traditional K for haplotypes
  const int S = eHapsCurrent_tc.n_slices;
  const int nSNPs = eHapsCurrent_tc.n_cols;
  //
  // new variables
  //
  if (snp_start_1_based == -1) {
      snp_start_1_based = 1;
      snp_end_1_based = nSNPs;
  }
  // variables working on nReads
  if (!pass_in_alphaBeta) {
      alphaHat_t = arma::zeros(K, nGrids);
      betaHat_t = arma::zeros(K, nGrids); 
      gamma_t = arma::zeros(K, nGrids);
      eMatGrid_t = arma::ones(K, nGrids);
  }
  arma::rowvec c = arma::zeros(1, nGrids);
  arma::mat eMatRead_t = arma::ones(K, nReads);
  arma::mat eMatHapOri_t;
  // for hapDosage, specifical definition of SNPs
  const int nSNPs_local = snp_end_1_based - snp_start_1_based + 1;
  arma::rowvec hapDosage = arma::zeros(1, nSNPs_local);
  arma::rowvec hapDosage_local = arma::zeros(1, nSNPs_local);  
  if (run_pseudo_haploid) {
      eMatHapOri_t = arma::zeros(K, nReads);
  }
  if (!update_in_place) {
      gammaSum0_tc = arma::zeros(K, nSNPs, S);
      gammaSum1_tc = arma::zeros(K, nSNPs, S);      
      alphaMatSum_tc = arma::zeros(K, nGrids - 1, S);      
      hapSum_tc = arma::zeros(K, nGrids, S);
      priorSum_m = arma::zeros(K, S);
  }
  // variables for transition matrix and initialization
  // int variables and such
  int iGrid, k, s;
  Rcpp::List to_return;
  Rcpp::List list_of_gammaK_t;  
  //  
  Rcpp::NumericVector alphaStart, betaEnd;
  arma::mat alphaHatBlocks_t, betaHatBlocks_t;
  Rcpp::List list_of_eMatRead_t;
  Rcpp::List list_of_hapDosage;
  arma::vec pRgivenH1(nReads);
  arma::vec pRgivenH2(nReads);
  double g_temp = 0;
  //
  // everything works on s here
  //
  for(s = 0; s < S; s++) {
      // always reset - always passed in
      alphaHat_t.fill(0);
      betaHat_t.fill(0);
      eMatRead_t.fill(1);
      eMatGrid_t.fill(1);
      //
      // so only need eMatGrid_t marginally
      //
      next_section="Make eMatGrid_t for special reference";
      prev=print_times(prev, suppressOutput, prev_section, next_section);
      prev_section=next_section;
      rcpp_ref_make_eMatGrid_t(eMatGrid_t, reference_haps, non_NA_cols, eHapsCurrent_tc, grid, reference_phred, s, iSample, maxEmissionMatrixDifference, rescale_eMatGrid_t, bound_eMatGrid_t);
      eMatRead_t = eMatGrid_t;
      //
      //
      next_section="Forward recursion";
      prev=print_times(prev, suppressOutput, prev_section, next_section);
      prev_section=next_section;
      //
      Rcpp_run_forward_haploid(alphaHat_t, c, eMatGrid_t, alphaMatCurrent_tc, transMatRate_tc_H, priorCurrent_m, s, alphaStart, run_fb_subset);
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
      next_section="Update things";
      prev=print_times(prev, suppressOutput, prev_section, next_section);
      prev_section=next_section;
      //
      hapSum_tc.slice(s) += gamma_t;
      priorSum_m.col(s) += gamma_t.col(0);
      //
      make_haploid_gammaUpdate_t(s, gammaSum0_tc, gammaSum1_tc, sampleReads, gamma_t, eHapsCurrent_tc, eMatRead_t, eMatHapOri_t, pRgivenH1, pRgivenH2, run_pseudo_haploid);
      //
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
