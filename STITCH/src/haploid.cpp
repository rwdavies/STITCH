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
    const int grid_offset = 0
);

Rcpp::List rcpp_make_fb_snp_offsets(
    const arma::mat& alphaHat_t,
    const arma::mat& betaHat_t,
    const arma::mat& blocks_for_output
);





//' @export
// [[Rcpp::export]]
void Rcpp_run_forward_haploid(
    arma::mat& alphaHat_t,
    arma::rowvec& c,
    const arma::mat& eMatHapSNP_t,
    const arma::cube& alphaMatCurrent_tc,
    const arma::cube& transMatRate_tc_H,
    const arma::mat& priorCurrent_m,
    const int s,
    const Rcpp::NumericVector alphaStart = 0,
    bool run_fb_subset = false
) {
    const int K = alphaMatCurrent_tc.n_rows;
    const int nGrids = alphaMatCurrent_tc.n_cols + 1;    
    //
    // initialize
    //
    int k;
    if (run_fb_subset == false) {
        for(k = 0; k < K; k++) {
            alphaHat_t(k, 0) = priorCurrent_m(k, s) * eMatHapSNP_t(k, 0);
        }
    } else {
        for(k=0; k < K; k++) {
            alphaHat_t(k, 0) = alphaStart(k);
        }
    }
    c(0) = 1 / sum(alphaHat_t.col(0));
    alphaHat_t.col(0) = alphaHat_t.col(0) * c(0);
    //
    //
    for(int iGrid = 1; iGrid < nGrids; iGrid++) {
        // NOTE - previously used this code, which is mathematically right
        // BUT since scaling is being used here, arma::sum(alphaHat_t.col(t - 1) is equal to 1 by definition
        // so can use the below (uncommented) code to simplify
        // alphaConst = transMatRate_t_H(1, t-1) * arma::sum(alphaHat_t.col(t - 1));
        //
        alphaHat_t.col(iGrid) = eMatHapSNP_t.col(iGrid) % (		   \
            transMatRate_tc_H(0, iGrid - 1, s) * alphaHat_t.col(iGrid - 1) + \
            transMatRate_tc_H(1, iGrid - 1, s) * alphaMatCurrent_tc.slice(s).col(iGrid - 1) );
        c(iGrid) = 1 / arma::sum(alphaHat_t.col(iGrid));
        alphaHat_t.col(iGrid) *= c(iGrid);
    }
    return ;
}


//' @export
// [[Rcpp::export]]
void Rcpp_run_backward_haploid(
    arma::mat& betaHat_t,
    arma::rowvec& c,
    const arma::mat& eMatHapSNP_t,
    const arma::cube& alphaMatCurrent_tc,
    const arma::cube& transMatRate_tc_H,
    const int s
) {
    const int nGrids = eMatHapSNP_t.n_cols;
    double x;
    arma::colvec e_times_b;
    for(int iGrid = nGrids - 2; iGrid >= 0; --iGrid) {
        e_times_b = eMatHapSNP_t.col(iGrid + 1) % betaHat_t.col(iGrid + 1);
        x = transMatRate_tc_H(1, iGrid, s) * sum(alphaMatCurrent_tc.slice(s).col(iGrid) % e_times_b);
        betaHat_t.col(iGrid) = c(iGrid) * (x + transMatRate_tc_H(0, iGrid, s) * e_times_b);
    }
    return;
}




//' @export
// [[Rcpp::export]]
void rcpp_make_eMatHap_t(
    arma::mat& eMatHap_t,
    const Rcpp::List& sampleReads,
    const arma::cube& eHapsCurrent_tc,
    const int s,
    const double maxDifferenceBetweenReads,
    const int Jmax,
    arma::mat& eMatHapOri_t,
    const arma::vec& pRgivenH1,
    const arma::vec& pRgivenH2,
    const bool run_pseudo_haploid = false,
    const bool rescale_eMatHap_t = true
) {
    //
    // constants
    //
    const int K = eHapsCurrent_tc.n_rows; // traditional K for haplotypes
    const int nReads = sampleReads.size();
    //
    // new variables
    //
    double pR, pA, eps, x;
    int j, k, J, readSNP, jj, iRead;
    //arma::mat eMatHap_t = arma::ones(K,nReads);
    //
    // now build
    //
    for(iRead = 0; iRead < nReads; iRead++) {
        Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
        // recal that below is what is used to set each element of sampleRead
        // note - this is no longer quite accurate
        // sampleReads.push_back(Rcpp::List::create(nU,d,phiU,pRU));
        J = as<int>(readData[0]); // number of Unique SNPs on read
        readSNP = as<int>(readData[1]); // leading SNP from read
        arma::ivec bqU = as<arma::ivec>(readData[2]); // bq for each SNP
        arma::ivec pRU = as<arma::ivec>(readData[3]); // position of each SNP from 0 to T-1
        // once each SNP is done, have P(read | k), can multiply to get P(read|(k1,k2))
        if(J >= Jmax) {
            J = Jmax;
        }
        for(j = 0; j <= J; j++) {
            if(bqU(j) < 0) {
                eps = pow(10,(double(bqU(j)) / 10));
                pR = 1 - eps;
                pA = eps / 3;
            }
            if(bqU(j) > 0) {
                eps = pow(10, (-double(bqU(j)) / 10));
                pR = eps / 3;
                pA = 1 - eps;
            }
            jj=pRU(j);
            eMatHap_t.col(iRead) %= ( eHapsCurrent_tc.slice(s).col(jj) * pA + (1 - eHapsCurrent_tc.slice(s).col(jj)) * pR);
            //
            if (run_pseudo_haploid == true) {
                x = pRgivenH1(iRead) / (pRgivenH1(iRead) + pRgivenH2(iRead));
                //
                for(k = 0; k < K; k++) {
                    // inefficient?
                    eMatHapOri_t(k, iRead) = eMatHap_t(k, iRead);
                    eMatHap_t(k,iRead) = x * eMatHap_t(k,iRead) + (1-x) * pRgivenH2(iRead);
                }
            }
        }
        //
        // cap P(read|k) to be within maxDifferenceBetweenReads orders of magnitude
        //
        if (rescale_eMatHap_t) {
            x=0;
            for(k=0; k<=K-1; k++)
                if(eMatHap_t(k,iRead)>x)
                    x=eMatHap_t(k,iRead);
            x = x / maxDifferenceBetweenReads;
            // x is the maximum now
            for(k=0; k<=K-1; k++)
                if(eMatHap_t(k,iRead)<x)
                    eMatHap_t(k,iRead) = x;
        }
    }
    return;
}
    


//' @export
// [[Rcpp::export]]
void rcpp_make_eMatHapSNP_t(
    arma::mat& eMatHapSNP_t,
    const arma::mat& eMatHap_t,
    const Rcpp::IntegerVector& H,
    const Rcpp::List sampleReads,
    const int hap,
    const int nGrids,
    const int run_fb_grid_offset = 0,
    const bool use_all_reads = false,
    const bool bound = false,
    const double maxEmissionMatrixDifference = 1000,
    const bool rescale = false
) {
    //
    int nReads = sampleReads.size(); //
    const int K = eMatHap_t.n_rows; // traditional K for haplotypes        
    // arma::mat eMatHapSNP_t = arma::ones(K, nGrids); // why is this called SNP? 
    int iRead, k, w, readSNP;
    bool proceed;
    //
    for(iRead = 0; iRead < nReads; iRead++) {
        proceed = false;
        if (use_all_reads) {
            proceed = true;
        } else {
            if (H(iRead) == hap) {
                proceed = true;
            }
        }
        if (proceed) {
            //w = wif(iRead) - run_fb_grid_offset;            
            // unclear how slow this is, but carry around less in RAM
            Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
            readSNP = as<int>(readData[1]); // leading SNP from read
            w = readSNP - run_fb_grid_offset;
            eMatHapSNP_t.col(w) %= eMatHap_t.col(iRead);
            //std::cout << "eMatHap_t.col(iRead) = " << eMatHap_t(0, iRead) << ", " << eMatHap_t(1, iRead) << ", " << eMatHap_t(2, iRead) << ", " << eMatHap_t(3, iRead);
            //                std::cout << std::endl;                
        }
    }
    // now - afterward - cap eMatHapSNP
    double x, rescale_val, d2;
    int t;
    if (bound) {
        for(t = 0; t < nGrids; t++) {
            // if this is less than exactly 1 (i.e. there are results here), proceed
            if (eMatHapSNP_t(0, t) < 1) {
                x = 0;
                for (k = 0; k < K; k++) {
                    if (eMatHapSNP_t(k, t) > x) {
                        x = eMatHapSNP_t(k, t);
                    }
                }
                // x is the maximum now
                rescale_val = 1 / x;        
                for (k = 0; k < K; k++) {
                    if (rescale) {                    
                        eMatHapSNP_t(k, t) *= rescale_val;
                    }
                    d2 = 1 / maxEmissionMatrixDifference;
                    if(eMatHapSNP_t(k, t) < (d2)) {
                        eMatHapSNP_t(k, t) = d2;
                    }
                }
            }
        }
    }
    return;
}


//' @export
// [[Rcpp::export]]
void make_haploid_gammaUpdate_t(
    int s,
    arma::cube& gammaUpdate_t,
    const Rcpp::List& sampleReads,
    const arma::mat& gamma_t,
    const arma::cube& eHapsCurrent_tc,
    const arma::mat& eMatHap_t,    
    const arma::mat& eMatHapOri_t,
    const arma::vec& pRgivenH1,
    const arma::vec& pRgivenH2,
    const bool run_pseudo_haploid = false
) {
    //
    //const int nSNPs = eHapsCurrent_t.n_cols;
    const int nReads = sampleReads.size();
    //
    //arma::cube gammaUpdate_t = arma::zeros(K, nSNPs, 2);
    int iRead, J, cr, t, j;
    Rcpp::List readData;
    arma::ivec bqU, pRU;
    arma::colvec gamma_t_col;
    //double d3, eps, a1, a2, y, d, d1, d2, b;
    double eps;
    double d3 = 1;
    double pR = -1;
    double pA = -1;
    int cr_prev = -1;
    arma::colvec a1_col, a2_col, y_col, d_col;
    arma::colvec b_col, d1_col, d2_col;
    //
    for(iRead = 0; iRead < nReads; iRead++) {
        readData = as<Rcpp::List>(sampleReads[iRead]);
        J = as<int>(readData[0]); // number of SNPs on read
        cr = as<int>(readData[1]); // central SNP in read
        bqU = as<arma::ivec>(readData[2]); // bq for each SNP
        pRU = as<arma::ivec>(readData[3]); // position of each SNP from 0 to T-1
        if (run_pseudo_haploid) {      
          d3 = pRgivenH1(iRead) / (pRgivenH1(iRead) + pRgivenH2(iRead));
        }
        if (cr > cr_prev) {
            // do not always need to update
            gamma_t_col = gamma_t.col(cr);
        }
        for(j = 0; j <= J; j++) {
            t=pRU(j);
            if(bqU(j)<0) {
                eps = pow(10,(double(bqU(j))/10));
                pR=1-eps;
                pA=eps/3;
            }
            if(bqU(j)>0) {
                eps = pow(10,(-double(bqU(j))/10));
                pR=eps/3;
                pA=1-eps;
            }
            //
            a1_col = pA * eHapsCurrent_tc.slice(s).col(t);
            a2_col = pR * (1 - eHapsCurrent_tc.slice(s).col(t));
            y_col = a1_col + a2_col;
            if (run_pseudo_haploid) {
                b_col = d3 * (eMatHapOri_t.col(iRead) / y_col);
                d1_col = a1_col % b_col;
                d2_col = a2_col % b_col;
                d_col = gamma_t_col / eMatHap_t.col(iRead);
                gammaUpdate_t.slice(0).col(t) += d_col % (d1_col);
                gammaUpdate_t.slice(0).col(t) += d_col % (d1_col + d2_col);
            } else {
                d_col = gamma_t_col / y_col;
                gammaUpdate_t.slice(0).col(t) += a1_col * d_col;
                gammaUpdate_t.slice(1).col(t) += y_col * d_col;
            }
        }
    }
    return;
}







//' @export
// [[Rcpp::export]]
Rcpp::List pseudoHaploid_update_model_9(const arma::vec& pRgivenH1, const arma::vec& pRgivenH2, const arma::mat& eMatHap_t1, const arma::mat& eMatHap_t2, const arma::mat& gamma_t1, const arma::mat& gamma_t2, const int K, const arma::ivec& srp) {
    // new stuff
    arma::vec pRgivenH1_new = arma::zeros(pRgivenH1.n_elem);
    arma::vec pRgivenH2_new = arma::zeros(pRgivenH2.n_elem);
    int k, t;
    double d1, d2, x1, x2;
    //
    for(std::size_t i_read=0; i_read < srp.n_elem; i_read++) {
        t = srp(i_read);
        x2 = pRgivenH2(i_read) / \
            (pRgivenH1(i_read) + pRgivenH2(i_read));
        x1 = 1 - x2;
        d1 = x2 * pRgivenH2(i_read);
        d2 = x1 * pRgivenH1(i_read);
        for(k=0; k < K; k++) {
            pRgivenH1_new(i_read) = pRgivenH1_new(i_read) +      \
                gamma_t1(k, t) * \
                (eMatHap_t1(k, i_read) - d1) / x1;
            pRgivenH2_new(i_read) = pRgivenH2_new(i_read) +  \
                gamma_t2(k, t) * \
                (eMatHap_t2(k, i_read) - d2) / x2;
        }
    }
    return List::create(
        Rcpp::Named("pRgivenH1_new") = pRgivenH1_new,
        Rcpp::Named("pRgivenH2_new") = pRgivenH2_new
                        );
}
    




void perform_haploid_per_sample_updates(
    int s,
    arma::mat& betaHat_t,
    arma::mat& eMatHapSNP_t,
    const arma::cube& transMatRate_tc_H,
    const arma::cube& alphaMatCurrent_tc,
    arma::cube& hapSum_tc,
    arma::mat& gamma_t,
    const bool update_in_place,
    const arma::mat& priorCurrent_m,
    arma::mat& priorSum_m,
    arma::cube& gammaUpdate_t,
    arma::cube& jUpdate_tc,
    const Rcpp::List& sampleReads,
    const arma::cube& eHapsCurrent_tc,
    const arma::mat& eMatHap_t,    
    const arma::mat& eMatHapOri_t,
    const arma::vec& pRgivenH1,
    const arma::vec& pRgivenH2,
    const bool run_pseudo_haploid,
    double prev,
    int suppressOutput,
    std::string prev_section,
    std::string next_section
) {
    //
    int k, iGrid;
    const int nGrids = alphaMatCurrent_tc.n_cols + 1;
    const int K = eHapsCurrent_tc.n_rows;
    //
    if (update_in_place) {  
        hapSum_tc.slice(s) += gamma_t;
        for(k = 0; k < K; k++) {
            priorSum_m(k, s) = priorSum_m(k, s) + gamma_t(k, 0);
        }
    }
    //
    // hap probs are gamma!
    //
    next_section="Gamma update";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //
    make_haploid_gammaUpdate_t(s, gammaUpdate_t, sampleReads, gamma_t, eHapsCurrent_tc, eMatHap_t, eMatHapOri_t, pRgivenH1, pRgivenH2, run_pseudo_haploid);
    //
    // make jUpdate
    //
    next_section="Make jUpdate";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //
    for(iGrid = 0; iGrid < nGrids - 1; iGrid++) {
        jUpdate_tc.slice(s).col(iGrid) += transMatRate_tc_H(1, iGrid, s) * (alphaMatCurrent_tc.slice(s).col(iGrid) % betaHat_t.col(iGrid + 1) % eMatHapSNP_t.col(iGrid + 1));
    }
}



//  calculate, multiply by 1/s, add to hapDosage
//  steal from rcpp_calculate_fbd_dosage(
void rcpp_calculate_hapDosage(
    arma::rowvec hapDosage
) {
    // do this
    // do something with this
    return;
}
    
//' @export
// [[Rcpp::export]]
Rcpp::List forwardBackwardHaploid(
    const Rcpp::List& sampleReads,
    const arma::cube& eHapsCurrent_tc,
    const arma::cube& alphaMatCurrent_tc,
    const arma::cube& transMatRate_tc_H,
    const arma::mat& priorCurrent_m,
    arma::mat& alphaHat_t,
    arma::mat& betaHat_t,
    arma::mat& gamma_t,    
    const double maxDifferenceBetweenReads,
    const double maxEmissionMatrixDifference,
    const int Jmax,
    const int suppressOutput,
    const int model,
    arma::cube& gammaUpdate_t,
    arma::cube& jUpdate_tc,
    arma::cube& hapSum_tc,
    arma::mat& priorSum_m,
    const arma::vec& pRgivenH1,
    const arma::vec& pRgivenH2,
    const bool run_pseudo_haploid,
    const arma::mat& blocks_for_output,
    const bool generate_fb_snp_offsets = false,
    const Rcpp::NumericVector alphaStart = 0,
    const Rcpp::NumericVector betaEnd = 0,
    const bool run_fb_subset = false,
    const int run_fb_grid_offset = 0, // this is 0-based
    const bool return_extra = false,
    const bool update_in_place = false,
    const bool pass_in_alphaBeta = false, // whether to pass in pre-made alphaHat, betaHat
    const bool output_haplotype_dosages = false, // whether to output state probabilities
    const int snp_start_1_based = -1,
    const int snp_end_1_based = -1,
    const Rcpp::IntegerVector grid = 0,
    const bool rescale = false
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
  // variables working on nReads
  if (!pass_in_alphaBeta) {
      alphaHat_t = arma::zeros(K, nGrids);
      betaHat_t = arma::zeros(K, nGrids); 
      gamma_t = arma::zeros(K, nGrids);
  }
  arma::rowvec c = arma::zeros(1, nGrids);
  arma::mat eMatHap_t = arma::ones(K, nReads);
  arma::mat eMatHapSNP_t = arma::ones(K, nGrids);
  arma::mat eMatHapOri_t;
  arma::rowvec hapDosage = arma::zeros(1, nSNPs);
  if (run_pseudo_haploid) {
      eMatHapOri_t = arma::zeros(K, nReads);
  }
  if (!update_in_place) {        
      gammaUpdate_t = arma::zeros(K, nSNPs, 2);
      jUpdate_tc = arma::zeros(K, nGrids - 1, S);
  }
  // variables for transition matrix and initialization
  // int variables and such
  int iGrid, k, s;
  Rcpp::List alphaBetaBlocks;
  Rcpp::List to_return;
  //
  // everything works on s here
  //
  for(s = 0; s < S; s++) {
      if (s > 0) {
          alphaHat_t.fill(0);
          betaHat_t.fill(0);
          eMatHap_t.fill(1);
          eMatHapSNP_t.fill(1);
      }
      //
      // eMatHap
      //
      next_section="Make eMatHap";
      prev=print_times(prev, suppressOutput, prev_section, next_section);
      prev_section=next_section;
      //
      // define eMatHap to work on the number of reads
      rcpp_make_eMatHap_t(eMatHap_t, sampleReads, eHapsCurrent_tc, s, maxDifferenceBetweenReads, Jmax, eMatHapOri_t, pRgivenH1, pRgivenH2, run_pseudo_haploid);
      //
      // once we have eMatHap, ie probabilities from reads, make eMatHapSNPsfrom this
      //
      next_section="Make eMat";
      prev=print_times(prev, suppressOutput, prev_section, next_section);
      prev_section=next_section;
      rcpp_make_eMatHapSNP_t(eMatHapSNP_t, eMatHap_t, 1, sampleReads, 1, nGrids, run_fb_grid_offset, true, true, maxEmissionMatrixDifference, rescale);
      //
      // forward recursion
      //
      next_section="Forward recursion";
      prev=print_times(prev, suppressOutput, prev_section, next_section);
      prev_section=next_section;
      //
      Rcpp_run_forward_haploid(alphaHat_t, c, eMatHapSNP_t, alphaMatCurrent_tc, transMatRate_tc_H, priorCurrent_m, s, alphaStart, run_fb_subset);
      //
      // backward recursion
      //
      next_section="Backward recursion";
      prev=print_times(prev, suppressOutput, prev_section, next_section);
      prev_section=next_section;
      if (run_fb_subset == false) {    
          betaHat_t.col(nGrids - 1).fill(c(nGrids - 1));
      } else {
          for(k = 0; k < K; k++) {
              betaHat_t(k, nGrids - 1) = betaEnd(k);
          }
      }
      Rcpp_run_backward_haploid(betaHat_t, c, eMatHapSNP_t, alphaMatCurrent_tc, transMatRate_tc_H, s);
      //
      // make gamma
      //
      next_section="Make gamma";
      prev=print_times(prev, suppressOutput, prev_section, next_section);
      prev_section=next_section;
      //
      gamma_t = alphaHat_t % betaHat_t;
      // normalize as well
      double g_temp = 0;
      for(iGrid = 0; iGrid < nGrids; iGrid++) {
          g_temp = 1 / c(iGrid);
          gamma_t.col(iGrid) *= g_temp;
      }
      //
      //
      if (output_haplotype_dosages) {
          // yuck - only with s = 1?
          arma::mat gammaEK_t = make_gammaEK_t_from_gammaK_t(
              gamma_t, K, grid,
              snp_start_1_based, snp_end_1_based, run_fb_grid_offset
          );
          to_return.push_back(gammaEK_t, "gammaEK_t");      
      }
      //
      // FIX THIS!
      //
      rcpp_calculate_hapDosage(hapDosage);
      //
      // if not fb subset, do updates
      //
      // three options
      // generate_fb_snp_offsets
      // run_fb_offset (do nothing - only hapDosage above)
      // normal (do updates)
      if (generate_fb_snp_offsets) {
          // will need to modify
          alphaBetaBlocks = rcpp_make_fb_snp_offsets(alphaHat_t, betaHat_t, blocks_for_output);
      } else if (!run_fb_subset) {
          //
          perform_haploid_per_sample_updates(s, betaHat_t, eMatHapSNP_t, alphaMatCurrent_tc, transMatRate_tc_H, hapSum_tc, gamma_t, update_in_place, priorCurrent_m, priorSum_m, gammaUpdate_t, jUpdate_tc, sampleReads, eHapsCurrent_tc, eMatHap_t, eMatHapOri_t, pRgivenH1, pRgivenH2, run_pseudo_haploid, prev, suppressOutput, prev_section, next_section);
      }
  }
  if (run_fb_subset == true) {
      to_return.push_back(hapDosage, "hapDosage");
      return(to_return);
  }
  //
  //
  //
  next_section="Done";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  //
  to_return.push_back(eMatHap_t, "eMatHap_t");
  to_return.push_back(eMatHapOri_t, "eMatHapOri_t");
  if (!update_in_place) {
      to_return.push_back(gammaUpdate_t, "gammaUpdate_t");
      to_return.push_back(jUpdate_tc, "jUpdate_tc");
  }
  if (return_extra) {
      to_return.push_back(alphaHat_t, "alphaHat_t");
      to_return.push_back(betaHat_t, "betaHat_t");
      to_return.push_back(eMatHapSNP_t, "eMatHapSNP_t");      
  }
  if (generate_fb_snp_offsets == true) {
      to_return.push_back(alphaBetaBlocks, "alphaBetaBlocks");
  }
  return(wrap(to_return));  
}


arma::ivec sample_path_from_alphaHat_t(const arma::mat & alphaHat_t, const arma::mat & transMatRate_t, const arma::mat & eMatHapSNP_t, const arma::mat & alphaMat_t, const int T, const int K, const arma::rowvec & c) {
    //
    arma::ivec sampled_state(T);
    int k, t;
    double sum = 0;
    double rand_uniform = 0;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    // somehow this makes dis(gen) give me random 0-1 doubles I think
    //
    // choose initial state
    //
    k = 0;  
    rand_uniform = dis(gen);
    sum = 0;
    while(k < K) {
        sum = sum + alphaHat_t(k, T - 1);
        if (sum > rand_uniform) {
            sampled_state(T - 1) = k;
            k = K;
        }
        k = k + 1;    
    }
    //
    // cycle
    //
    double samp_vector_sum;
    double norm;
    arma::rowvec samp_vector = arma::zeros(1, K);  
    int prev_state;
    //for(t = T - 2; t >= 0; t--) {
    for(t = T - 2; t >= 0; t--) {
        prev_state = sampled_state(t + 1);
        samp_vector.fill(transMatRate_t(1, t) * alphaMat_t(prev_state, t));
        samp_vector(prev_state)=samp_vector(prev_state) + transMatRate_t(0, t);
        samp_vector_sum = 0;
        for(k=0; k<=K-1; k++) {
            samp_vector(k)=samp_vector(k) * alphaHat_t(k, t);
            samp_vector_sum=samp_vector_sum + samp_vector(k);
        }
        norm = alphaHat_t(prev_state, t + 1) / eMatHapSNP_t(prev_state, t + 1) / c(t + 1);
        // if ((pow(samp_vector_sum / norm, 2) - 1) > 1e-4) {
        //     std::cout << "BAD COUNT on t=" << t << std::endl;
        //     std::cout << "samp_vector_sum=" << samp_vector_sum << std::endl;
        //     std::cout << "norm=" << norm << std::endl;
        //     return(sampled_state);          
        // }
        // sample state using crazy C++
        k = 0;  
        rand_uniform = dis(gen);
        sum = 0;
        while(k < K) {
            sum = sum + samp_vector(k) / norm;
            if (sum > rand_uniform) {
                sampled_state(t) = k;
                k = K;
            }
            k = k + 1;    
        }
    }
    return(sampled_state);
}



//' @export
// [[Rcpp::export]]
arma::ivec rcpp_sample_path(const arma::rowvec read_labels, const arma::mat eMatHap_t, const Rcpp::List& sampleReads, const int nReads, const arma::mat& eHaps_t, const double maxDifferenceBetweenReads, const int Jmax, const arma::vec pi, const arma::mat & transMatRate_t_H, const arma::mat& alphaMat_t) {
  //
  // constants
  //
  const int T = eHaps_t.n_cols;
  const int K = eHaps_t.n_rows; // traditional K for haplotypes
  //
  // new variables
  //
  int iRead, readSNP, k, k1, t;
  double alphaConst;
  arma::mat eMatHapSNP_t = arma::ones(K,T);
  arma::mat alphaHat_t = arma::zeros(K,T);
  arma::rowvec c = arma::zeros(1,T);
  //
  //
  // initialize emission matrix from eMatHapSNP_t which is constant
  //
  //
  for(iRead=0; iRead<=nReads-1; iRead++) {
      if (read_labels(iRead) == 1) {
          Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
          readSNP = as<int>(readData[1]); // leading SNP from read
          for(k=0; k<=K-1; k++) {
              eMatHapSNP_t(k, readSNP) = eMatHapSNP_t(k, readSNP) *     \
                  eMatHap_t(k,iRead);
          }
      }
  }
  //
  //
  // forward algorithm
  //
  //
  for(k1=0; k1<=K-1; k1++)
      alphaHat_t(k1,0) = pi(k1) * eMatHapSNP_t(k1,0);
  c(0) = 1 / sum(alphaHat_t.col(0));
  alphaHat_t.col(0) = alphaHat_t.col(0) * c(0);
  //
  for(t=1; t<=T-1; t++) {
      // make the constant
      alphaConst=0;
      for(k=0; k<=K-1; k++)
          alphaConst = alphaConst + alphaHat_t(k,t-1);
      alphaConst = alphaConst * transMatRate_t_H(1, t-1);
      //
      // each entry is emission * (no change * that value + constant)
      for(k=0; k<=K-1; k++)
          alphaHat_t(k,t) =  eMatHapSNP_t(k,t) *  \
              ( transMatRate_t_H(0, t-1) * alphaHat_t(k,t-1) +   \
                alphaConst * alphaMat_t(k,t-1));
      //
      c(t) = 1 / sum(alphaHat_t.col(t));
      alphaHat_t.col(t) = alphaHat_t.col(t) * c(t);
  }
  //
  //
  // sample a path
  //
  //
  arma::ivec to_out = sample_path_from_alphaHat_t(alphaHat_t, transMatRate_t_H, eMatHapSNP_t, alphaMat_t, T, K, c);
  return(to_out);
}



      




//' @export
// [[Rcpp::export]]
arma::mat rcpp_calculate_many_likelihoods(const arma::mat swap_mat, const Rcpp::List reads_at_SNPs, const arma::mat eMatHap_t, const Rcpp::List& sampleReads, const int nReads, const arma::mat& eHaps_t, const double maxDifferenceBetweenReads, const int Jmax, const arma::vec pi, const arma::mat& transMatRate_t, const arma::mat& alphaMat_t) {
  //
  // constants
  //
  const int nSwap = swap_mat.n_rows;
  const int T = eHaps_t.n_cols;
  const int K = eHaps_t.n_rows; // traditional K for haplotypes
  //
  // new variables
  //
  int iRead, readSNP, k, t, s, i;
  double alphaConst_hap1, alphaConst_hap2, c_hap1, c_hap2;
  arma::mat eMatHapSNPSwap_hap1 = arma::ones(nSwap, K);
  arma::mat eMatHapSNPSwap_hap2 = arma::ones(nSwap, K);  
  arma::mat alphaHatSwapPresent_hap1 = arma::zeros(nSwap, K);
  arma::mat alphaHatSwapPresent_hap2 = arma::zeros(nSwap, K);  
  arma::mat alphaHatSwapFuture_hap1 = arma::zeros(nSwap, K);
  arma::mat alphaHatSwapFuture_hap2 = arma::zeros(nSwap, K);    
  arma::mat cFinal = arma::zeros(nSwap, 2);
  //
  //
  //
  //
  //
  Rcpp::NumericVector reads_at_SNP = as<Rcpp::NumericVector>(reads_at_SNPs[0]);
  eMatHapSNPSwap_hap1.fill(1);
  eMatHapSNPSwap_hap2.fill(1);    
  if (reads_at_SNP(0) >= 0) {
      for(i = 0; i < reads_at_SNP.size(); i++) {
          iRead = reads_at_SNP(i);
          Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
          readSNP = as<int>(readData[1]); // leading SNP from read
          for(s=0; s <= nSwap - 1; s++) {
              if (swap_mat(s, iRead) == 1) {
                  for(k=0; k<=K-1; k++) {                  
                      eMatHapSNPSwap_hap1(s, k) = eMatHapSNPSwap_hap1(s, k) * eMatHap_t(k, iRead);
                  }
              } else {
                  for(k=0; k<=K-1; k++) {                  
                      eMatHapSNPSwap_hap2(s, k) = eMatHapSNPSwap_hap2(s, k) * eMatHap_t(k, iRead);
                  }
              }
          }
      }
  }
  //
  // rest of initialization
  //
  for(s=0; s <= nSwap - 1; s++) {
      for(k=0; k <= K-1; k++) {
          alphaHatSwapFuture_hap1(s, k) = pi(k) * eMatHapSNPSwap_hap1(s, k);
          alphaHatSwapFuture_hap2(s, k) = pi(k) * eMatHapSNPSwap_hap2(s, k);          
      }
      c_hap1 = 1 / sum(alphaHatSwapFuture_hap1.row(s));
      c_hap2 = 1 / sum(alphaHatSwapFuture_hap2.row(s));      
      alphaHatSwapFuture_hap1.row(s) = alphaHatSwapFuture_hap1.row(s) * c_hap1;
      alphaHatSwapFuture_hap2.row(s) = alphaHatSwapFuture_hap2.row(s) * c_hap2;      
      cFinal(s, 0) = cFinal(s, 0) + log(c_hap1);
      cFinal(s, 1) = cFinal(s, 1) + log(c_hap2);      
  }
  //
  //
  //
  for(t=1; t<=T-1; t++) {
      //std::cout << "SNP t= " << t << "\n";          
      // re-set Present one, can ignore future, is overridden
      alphaHatSwapPresent_hap1 = alphaHatSwapFuture_hap1;
      alphaHatSwapPresent_hap2 = alphaHatSwapFuture_hap2;      
      //
      // calculate effect of reads, if relevant
      //
      //std::cout << "read part\n";        
      Rcpp::NumericVector reads_at_SNP = as<Rcpp::NumericVector>(reads_at_SNPs[t]);
      eMatHapSNPSwap_hap1.fill(1);
      eMatHapSNPSwap_hap2.fill(1);      
      if (reads_at_SNP(0) >= 0) {
          for(i = 0; i < reads_at_SNP.size(); i++) {
              iRead = reads_at_SNP(i);
              Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
              readSNP = as<int>(readData[1]); // leading SNP from read
              for(s=0; s <= nSwap - 1; s++) {
                  if (swap_mat(s, iRead) == 1) {                  
                      for(k=0; k<=K-1; k++) {
                          eMatHapSNPSwap_hap1(s, k) = eMatHapSNPSwap_hap1(s, k) * eMatHap_t(k, iRead);
                      }
                  } else {
                      for(k=0; k<=K-1; k++) {
                          eMatHapSNPSwap_hap2(s, k) = eMatHapSNPSwap_hap2(s, k) * eMatHap_t(k, iRead);
                      }
                  }
              }
          }
      }
      //
      // now apply
      //
      //std::cout << "apply\n";              
      for(s=0; s<=nSwap - 1; s++) {
          alphaConst_hap1 = sum(alphaHatSwapPresent_hap1.row(s)) * transMatRate_t(1, t-1);
          alphaConst_hap2 = sum(alphaHatSwapPresent_hap2.row(s)) * transMatRate_t(1, t-1);
          for(k=0; k<=K-1; k++) {
              alphaHatSwapFuture_hap1(s, k) = eMatHapSNPSwap_hap1(s, k) * \
                  ( transMatRate_t(0, t-1) * alphaHatSwapPresent_hap1(s, k) + \
                    alphaConst_hap1 * alphaMat_t(k,t-1));
              alphaHatSwapFuture_hap2(s, k) = eMatHapSNPSwap_hap2(s, k) * \
                  ( transMatRate_t(0, t-1) * alphaHatSwapPresent_hap2(s, k) + \
                    alphaConst_hap2 * alphaMat_t(k,t-1));
          }
          c_hap1 = 1 / sum(alphaHatSwapFuture_hap1.row(s));
          c_hap2 = 1 / sum(alphaHatSwapFuture_hap2.row(s));      
          alphaHatSwapFuture_hap1.row(s) = alphaHatSwapFuture_hap1.row(s) * c_hap1;
          alphaHatSwapFuture_hap2.row(s) = alphaHatSwapFuture_hap2.row(s) * c_hap2;      
          cFinal(s, 0) = cFinal(s, 0) + log(c_hap1);
          cFinal(s, 1) = cFinal(s, 1) + log(c_hap2);      
          //std::cout << "s=" << s << "\n";
          //std::cout << "c_hap1=" << c_hap1 << "\n";
          //std::cout << "c_hap2=" << c_hap2 << "\n";
      }
  }
  //
  //
  return(cFinal);
}

