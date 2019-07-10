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
    const arma::mat& eMatGrid_t,
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
            alphaHat_t(k, 0) = priorCurrent_m(k, s) * eMatGrid_t(k, 0);
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
        alphaHat_t.col(iGrid) = eMatGrid_t.col(iGrid) % (		   \
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
    const arma::mat& eMatGrid_t,
    const arma::cube& alphaMatCurrent_tc,
    const arma::cube& transMatRate_tc_H,
    const int s
) {
    const int nGrids = eMatGrid_t.n_cols;
    double x;
    arma::colvec e_times_b;
    for(int iGrid = nGrids - 2; iGrid >= 0; --iGrid) {
        e_times_b = eMatGrid_t.col(iGrid + 1) % betaHat_t.col(iGrid + 1);
        x = transMatRate_tc_H(1, iGrid, s) * sum(alphaMatCurrent_tc.slice(s).col(iGrid) % e_times_b);
        betaHat_t.col(iGrid) = c(iGrid) * (x + transMatRate_tc_H(0, iGrid, s) * e_times_b);
    }
    return;
}




//' @export
// [[Rcpp::export]]
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
) {
    next_section="make eMatRead_t";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
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
    //arma::mat eMatRead_t = arma::ones(K,nReads);
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
            eMatRead_t.col(iRead) %= ( eHapsCurrent_tc.slice(s).col(jj) * pA + (1 - eHapsCurrent_tc.slice(s).col(jj)) * pR);
            //
            if (run_pseudo_haploid == true) {
                x = pRgivenH1(iRead) / (pRgivenH1(iRead) + pRgivenH2(iRead));
                //
                for(k = 0; k < K; k++) {
                    // inefficient?
                    eMatHapOri_t(k, iRead) = eMatRead_t(k, iRead);
                    eMatRead_t(k,iRead) = x * eMatRead_t(k,iRead) + (1-x) * pRgivenH2(iRead);
                }
            }
        }
        //
        // cap P(read|k) to be within maxDifferenceBetweenReads orders of magnitude
        //
        if (rescale_eMatRead_t) {
            x=0;
            for(k=0; k<=K-1; k++)
                if(eMatRead_t(k,iRead)>x)
                    x=eMatRead_t(k,iRead);
            x = x / maxDifferenceBetweenReads;
            // x is the maximum now
            for(k=0; k<=K-1; k++)
                if(eMatRead_t(k,iRead)<x)
                    eMatRead_t(k,iRead) = x;
        }
    }
    return;
}
    


//' @export
// [[Rcpp::export]]
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
) {
    //
    next_section="Make eMat";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //    
    int nReads = sampleReads.size(); //
    const int K = eMatRead_t.n_rows; // traditional K for haplotypes        
    // arma::mat eMatGrid_t = arma::ones(K, nGrids); // why is this called SNP? 
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
            eMatGrid_t.col(w) %= eMatRead_t.col(iRead);
            //std::cout << "eMatRead_t.col(iRead) = " << eMatRead_t(0, iRead) << ", " << eMatRead_t(1, iRead) << ", " << eMatRead_t(2, iRead) << ", " << eMatRead_t(3, iRead);
            //                std::cout << std::endl;                
        }
    }
    // now - afterward - cap eMatHapSNP
    double x, rescale_val, d2;
    int t;
    if (bound) {
        for(t = 0; t < nGrids; t++) {
            // if this is less than exactly 1 (i.e. there are results here), proceed
            if (eMatGrid_t(0, t) < 1) {
                x = 0;
                for (k = 0; k < K; k++) {
                    if (eMatGrid_t(k, t) > x) {
                        x = eMatGrid_t(k, t);
                    }
                }
                // x is the maximum now
                rescale_val = 1 / x;        
                for (k = 0; k < K; k++) {
                    if (rescale) {                    
                        eMatGrid_t(k, t) *= rescale_val;
                    }
                    d2 = 1 / maxEmissionMatrixDifference;
                    if(eMatGrid_t(k, t) < (d2)) {
                        eMatGrid_t(k, t) = d2;
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
                d_col = gamma_t_col / eMatRead_t.col(iRead);
                gammaSum0_tc.slice(s).col(t) += d_col % (d1_col);
                gammaSum1_tc.slice(s).col(t) += d_col % (d1_col + d2_col);
            } else {
                d_col = gamma_t_col / y_col;
                gammaSum0_tc.slice(s).col(t) += a1_col % d_col;
                gammaSum1_tc.slice(s).col(t) += y_col % d_col;
            }
        }
    }
    return;
}







    




void perform_haploid_per_sample_updates(
    int s,
    arma::mat& betaHat_t,
    arma::mat& eMatGrid_t,
    const Rcpp::List& sampleReads,    
    const arma::cube& eHapsCurrent_tc,
    const arma::cube& alphaMatCurrent_tc,
    const arma::cube& transMatRate_tc_H,
    const arma::mat& priorCurrent_m,    
    const arma::mat& eMatRead_t,    
    const arma::mat& eMatHapOri_t,
    arma::mat& gamma_t,
    const bool update_in_place,
    arma::cube& gammaSum0_tc,
    arma::cube& gammaSum1_tc,
    arma::cube& alphaMatSum_tc,
    arma::mat& priorSum_m,
    const arma::vec& pRgivenH1,
    const arma::vec& pRgivenH2,
    const bool run_pseudo_haploid,
    double& prev,
    int suppressOutput,
    std::string& prev_section,
    std::string& next_section
) {
    //
    const int nGrids = alphaMatCurrent_tc.n_cols + 1;
    //
    priorSum_m.col(s) += gamma_t.col(0);
    //
    next_section="Gamma update";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //
    make_haploid_gammaUpdate_t(s, gammaSum0_tc, gammaSum1_tc, sampleReads, gamma_t, eHapsCurrent_tc, eMatRead_t, eMatHapOri_t, pRgivenH1, pRgivenH2, run_pseudo_haploid);
    //
    // make jUpdate
    //
    next_section="Make jUpdate (alphaMatSum_tc)";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //
    for(int iGrid = 0; iGrid < nGrids - 1; iGrid++) {
        alphaMatSum_tc.slice(s).col(iGrid) += transMatRate_tc_H(1, iGrid, s) * (alphaMatCurrent_tc.slice(s).col(iGrid) % betaHat_t.col(iGrid + 1) % eMatGrid_t.col(iGrid + 1));
        //       jUpdate_t.col(t) += transMatRate_t_H(1, t) * (alphaMat_t.col(t) % betaHat_t.col(t + 1) % eMatGrid_t.col(t + 1));
    }
    return;
}



//  calculate, multiply by 1/s, add to hapDosage
//  steal from rcpp_calculate_fbd_dosage(

//' @export
// [[Rcpp::export]]
arma::rowvec rcpp_calculate_hapDosage(
    const arma::cube& eHapsCurrent_tc,
    const int s,
    const arma::mat& gamma_t,
    const Rcpp::IntegerVector& grid,
    const int snp_start_1_based,
    const int snp_end_1_based,
    const int run_fb_grid_offset = 0
) {
    //
    const int nSNPs = snp_end_1_based - snp_start_1_based + 1;
    arma::rowvec hapDosage = arma::zeros(1, nSNPs);
    int iSNP, t, cur_grid;
    arma::colvec gamma_t_col;
    int prev_grid = -1;
    for(iSNP = 0; iSNP < nSNPs; iSNP++) {
        t = iSNP + snp_start_1_based - 1;
        cur_grid = grid(t) - run_fb_grid_offset;
        if (cur_grid > prev_grid) {
            gamma_t_col = gamma_t.col(cur_grid);
            prev_grid = cur_grid;
        }
        hapDosage(iSNP) += arma::sum(eHapsCurrent_tc.slice(s).col(t) % gamma_t_col);
    }
    return(hapDosage);
}





//' @export
// [[Rcpp::export]]
Rcpp::List pseudoHaploid_update_model_9(
    const arma::mat& pRgivenH1_m,
    const arma::mat& pRgivenH2_m,
    const Rcpp::List& list_of_eMatRead_t1,
    const Rcpp::List& list_of_eMatRead_t2,
    const Rcpp::List& list_of_gamma_t1,
    const Rcpp::List& list_of_gamma_t2,
    const int K,
    const arma::ivec& srp
) {
    // new stuff
    const int S = pRgivenH1_m.n_cols;
    //
    arma::mat pRgivenH1_m_new = arma::zeros(pRgivenH1_m.n_rows, S);
    arma::mat pRgivenH2_m_new = arma::zeros(pRgivenH2_m.n_rows, S);
    int k, t, s;
    double d1, d2, x1, x2;
    arma::mat gamma_t1, gamma_t2, eMatRead_t1, eMatRead_t2;
    //
    for(s = 0; s < S; s++) {
        //
        gamma_t1 = as<arma::mat>(list_of_gamma_t1[s]);
        gamma_t2 = as<arma::mat>(list_of_gamma_t2[s]);
        eMatRead_t1 = as<arma::mat>(list_of_eMatRead_t1[s]);
        eMatRead_t2 = as<arma::mat>(list_of_eMatRead_t2[s]);        
        //
        for(std::size_t i_read=0; i_read < srp.n_elem; i_read++) {
            t = srp(i_read);
            x2 = pRgivenH2_m(i_read, s) / \
                 (pRgivenH1_m(i_read, s) + pRgivenH2_m(i_read, s));
            x1 = 1 - x2;
            d1 = x2 * pRgivenH2_m(i_read, s);
            d2 = x1 * pRgivenH1_m(i_read, s);
            for(k=0; k < K; k++) {
                pRgivenH1_m_new(i_read, s) += gamma_t1(k, t) * (eMatRead_t1(k, i_read) - d1) / x1;
                pRgivenH2_m_new(i_read, s) += gamma_t2(k, t) * (eMatRead_t2(k, i_read) - d2) / x2;
            }
        }
    }
    return List::create(
        Rcpp::Named("pRgivenH1_m_new") = pRgivenH1_m_new,
        Rcpp::Named("pRgivenH2_m_new") = pRgivenH2_m_new
    );
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
    const arma::mat& pRgivenH1_m,
    const arma::mat& pRgivenH2_m,
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
  Rcpp::List list_of_gamma_t;
  //  
  Rcpp::NumericVector alphaStart, betaEnd;
  arma::mat alphaHatBlocks_t, betaHatBlocks_t;
  Rcpp::List alphaBetaBlocks;
  Rcpp::List list_of_alphaBetaBlocks;
  Rcpp::List list_of_eMatRead_t;
  Rcpp::List list_of_hapDosage;
  arma::vec pRgivenH1(nReads);
  arma::vec pRgivenH2(nReads);
  //
  // everything works on s here
  //
  for(s = 0; s < S; s++) {
      // reset if > first s, or if they're passed in
      if ((s > 0) | pass_in_alphaBeta) {
          alphaHat_t.fill(0);
          betaHat_t.fill(0);
          eMatRead_t.fill(1);
          eMatGrid_t.fill(1);
          if (run_pseudo_haploid) {
              eMatHapOri_t.fill(0);
          }
      }
      if (run_fb_subset) {
          alphaBetaBlocks = as<Rcpp::List>(prev_list_of_alphaBetaBlocks[s]);
          alphaHatBlocks_t = as<arma::mat>(alphaBetaBlocks["alphaHatBlocks_t"]);
          betaHatBlocks_t = as<arma::mat>(alphaBetaBlocks["betaHatBlocks_t"]);
          alphaStart = alphaHatBlocks_t.col(i_snp_block_for_alpha_beta);
          betaEnd = betaHatBlocks_t.col(i_snp_block_for_alpha_beta);
      }
      if (run_pseudo_haploid) {
          for(int iRead = 0; iRead < nReads; iRead++) {
              pRgivenH1(iRead) = pRgivenH1_m(iRead, s);
              pRgivenH2(iRead) = pRgivenH2_m(iRead, s);
          }
      }
      //
      rcpp_make_eMatRead_t(eMatRead_t, sampleReads, eHapsCurrent_tc, s, maxDifferenceBetweenReads, Jmax, eMatHapOri_t, pRgivenH1, pRgivenH2, prev, suppressOutput, prev_section, next_section, run_pseudo_haploid);
      //
      rcpp_make_eMatGrid_t(eMatGrid_t, eMatRead_t, 1, sampleReads, 1, nGrids, prev, suppressOutput, prev_section,next_section, run_fb_grid_offset, true, true, maxEmissionMatrixDifference, rescale);
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
      if (run_fb_subset == false) {    
          betaHat_t.col(nGrids - 1).fill(c(nGrids - 1));
      } else {
          for(k = 0; k < K; k++) {
              betaHat_t(k, nGrids - 1) = betaEnd(k);
          }
      }
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
      double g_temp = 0;
      for(iGrid = 0; iGrid < nGrids; iGrid++) {
          g_temp = 1 / c(iGrid);
          gamma_t.col(iGrid) *= g_temp;
      }
      //
      if (return_gamma | run_pseudo_haploid) {
          list_of_gamma_t.push_back(gamma_t, "gamma_t");      
      }
      //
      //
      if (output_haplotype_dosages) {
          // yuck - only with s = 1?
          arma::mat gammaEK_t = make_gammaEK_t_from_gammaK_t(
              gamma_t, K, grid,
              snp_start_1_based, snp_end_1_based,
              prev, suppressOutput, prev_section, next_section,
              run_fb_grid_offset
          );
          to_return.push_back(gammaEK_t, "gammaEK_t");      
      }
      //
      next_section="Make hap dosage";
      prev=print_times(prev, suppressOutput, prev_section, next_section);
      prev_section=next_section;
      //
      hapDosage_local = rcpp_calculate_hapDosage(eHapsCurrent_tc, s, gamma_t, grid, snp_start_1_based, snp_end_1_based, run_fb_grid_offset);
      hapDosage = hapDosage + hapDosage_local;
      if (return_hapDosage) {
          list_of_hapDosage.push_back(hapDosage_local);
      }
      //
      // if not fb subset, do updates
      //
      // three options
      // generate_fb_snp_offsets
      // run_fb_offset (do nothing - only hapDosage above)
      // normal (do updates)
      if (generate_fb_snp_offsets) {
          // will need to modify
          list_of_alphaBetaBlocks.push_back(rcpp_make_fb_snp_offsets(alphaHat_t, betaHat_t, blocks_for_output), "alphaBetaBlocks");
      }
      if (!run_fb_subset) {
          // hapSum here
          hapSum_tc.slice(s) += gamma_t;          
          if (!generate_fb_snp_offsets) {
              //
              perform_haploid_per_sample_updates(
                  s, betaHat_t, eMatGrid_t, sampleReads,    
                  eHapsCurrent_tc,  alphaMatCurrent_tc,  transMatRate_tc_H,  priorCurrent_m,    
                  eMatRead_t, eMatHapOri_t, gamma_t, update_in_place,
                  gammaSum0_tc, gammaSum1_tc, alphaMatSum_tc, priorSum_m,
                  pRgivenH1, pRgivenH2, run_pseudo_haploid,
                  prev, suppressOutput, prev_section, next_section
              );
          }
      }
      if (return_gammaK) {
          if (s == (S - 1)) {
              to_return.push_back(gamma_t, "gammaK_t"); // equivalent
          }
      }
      if (run_pseudo_haploid & !run_fb_subset & !generate_fb_snp_offsets) {
          list_of_eMatRead_t.push_back(eMatRead_t, "eMatRead_t");
      }
  }
  next_section="Done";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  if (return_gamma | run_pseudo_haploid) {
      to_return.push_back(list_of_gamma_t, "list_of_gamma_t");
  }
  if (run_pseudo_haploid) {
      to_return.push_back(eMatRead_t, "eMatRead_t");
      to_return.push_back(list_of_eMatRead_t , "list_of_eMatRead_t");
  }
  if (return_hapDosage) {
      if (S > 1) {
          // divide by S, here always just add
          hapDosage *= 1 / double(S);
      }
      to_return.push_back(hapDosage, "hapDosage");
      to_return.push_back(list_of_hapDosage, "list_of_hapDosage");
  }
  //
  if (generate_fb_snp_offsets) {
      to_return.push_back(list_of_alphaBetaBlocks, "list_of_alphaBetaBlocks");
  }
  if (run_fb_subset | generate_fb_snp_offsets) {
      return(to_return);
  }
  //
  if (!update_in_place) {
      to_return.push_back(gammaSum0_tc, "gammaSum0_tc");      
      to_return.push_back(gammaSum1_tc, "gammaSum1_tc");
      to_return.push_back(alphaMatCurrent_tc, "alphaMatCurrent_tc");
      to_return.push_back(hapSum_tc, "hapSum_tc");
      to_return.push_back(priorSum_m, "priorSum_m");      
  }
  if (return_extra) {
      to_return.push_back(alphaHat_t, "alphaHat_t");
      to_return.push_back(betaHat_t, "betaHat_t");
      to_return.push_back(eMatGrid_t, "eMatGrid_t");
      to_return.push_back(eMatRead_t, "eMatRead_t");
      to_return.push_back(eMatHapOri_t, "eMatHapOri_t");
  }
  next_section="Done";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  return(to_return);  
}


arma::ivec sample_path_from_alphaHat_t(const arma::mat & alphaHat_t, const arma::cube & transMatRate_tc, const arma::mat & eMatGrid_t, const arma::cube & alphaMatCurrent_tc, const int T, const int K, const arma::rowvec & c, const int s) {
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
        samp_vector.fill(transMatRate_tc(1, t, s) * alphaMatCurrent_tc(prev_state, t, s));
        samp_vector(prev_state)=samp_vector(prev_state) + transMatRate_tc(0, t, s);
        samp_vector_sum = 0;
        for(k=0; k<=K-1; k++) {
            samp_vector(k)=samp_vector(k) * alphaHat_t(k, t);
            samp_vector_sum=samp_vector_sum + samp_vector(k);
        }
        norm = alphaHat_t(prev_state, t + 1) / eMatGrid_t(prev_state, t + 1) / c(t + 1);
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
arma::ivec rcpp_sample_path(
    const arma::rowvec read_labels,
    const arma::mat eMatRead_t,
    const Rcpp::List& sampleReads,
    const double maxDifferenceBetweenReads,
    const int Jmax,
    const arma::mat& priorCurrent_m,
    const arma::cube& transMatRate_tc_H,
    const arma::cube& alphaMatCurrent_tc,
    const int s
) {
    //
    // cponstants
    //
    const int nGrids = alphaMatCurrent_tc.n_cols + 1;
    const int nReads = sampleReads.size();
    const int K = alphaMatCurrent_tc.n_rows; // traditional K for haplotypes
    //
    // new variables
    //
    int iRead, readSNP, k, k1, t;
    double alphaConst;
    arma::mat eMatGrid_t = arma::ones(K, nGrids);
    arma::mat alphaHat_t = arma::zeros(K, nGrids);
    arma::rowvec c = arma::zeros(1, nGrids);
    //
    //
    // initialize emission matrix from eMatGrid_t which is constant
    //
    //
    for(iRead=0; iRead<=nReads-1; iRead++) {
        if (read_labels(iRead) == 1) {
            Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
            readSNP = as<int>(readData[1]); // leading SNP from read
            for(k=0; k<=K-1; k++) {
                eMatGrid_t(k, readSNP) = eMatGrid_t(k, readSNP) * eMatRead_t(k,iRead);
            }
        }
    }
    //
    //
    // forward algorithm
    //
    //
    for(k1=0; k1<=K-1; k1++)
        alphaHat_t(k1,0) = priorCurrent_m(k1, s) * eMatGrid_t(k1,0);
    c(0) = 1 / sum(alphaHat_t.col(0));
    alphaHat_t.col(0) = alphaHat_t.col(0) * c(0);
    //
    for(t=1; t<= nGrids-1; t++) {
        // make the constant
        alphaConst=0;
        for(k=0; k<=K-1; k++)
            alphaConst = alphaConst + alphaHat_t(k,t-1);
        alphaConst = alphaConst * transMatRate_tc_H(1, t-1, s);
        //
        // each entry is emission * (no change * that value + constant)
        for(k=0; k<=K-1; k++)
            alphaHat_t(k,t) =  eMatGrid_t(k,t) *                 \
                ( transMatRate_tc_H(0, t-1, s) * alphaHat_t(k,t-1) +      \
                  alphaConst * alphaMatCurrent_tc(k,t-1, s));
        //
        c(t) = 1 / sum(alphaHat_t.col(t));
        alphaHat_t.col(t) = alphaHat_t.col(t) * c(t);
    }
    //
    //
    // sample a path
    //
    //
    arma::ivec to_out = sample_path_from_alphaHat_t(alphaHat_t, transMatRate_tc_H, eMatGrid_t, alphaMatCurrent_tc, nGrids, K, c, s);
    return(to_out);
}


      




