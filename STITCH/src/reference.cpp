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
    const int model,
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
    const bool rescale_eMatGrid_t = false
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
      rcpp_make_eMatRead_t(eMatRead_t, sampleReads, eHapsCurrent_tc, s, maxDifferenceBetweenReads, Jmax, eMatHapOri_t, pRgivenH1, pRgivenH2, prev, suppressOutput, prev_section, next_section, run_pseudo_haploid);
      //
      rcpp_make_eMatGrid_t(eMatGrid_t, eMatRead_t, 1, sampleReads, 1, nGrids, prev, suppressOutput, prev_section,next_section, run_fb_grid_offset, true, true, maxEmissionMatrixDifference, rescale_eMatGrid_t);
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
  //
  //
  next_section="Done";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  return(to_return);  
}
