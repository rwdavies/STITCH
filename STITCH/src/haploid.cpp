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
    const arma::mat& alphaMat_t,    
    const arma::mat& transMatRate_t_H,
    const int& T,
    const int& K,
    const arma::vec pi,
    const Rcpp::NumericVector alphaStart = 0,
    bool run_fb_subset = false
) {
    double alphaConst;
    //
    // initialize
    //
    int k;
    if (run_fb_subset == false) {
        for(k = 0; k < K; k++) {
            alphaHat_t(k,0) = pi(k) * eMatHapSNP_t(k,0);
        }
    } else {
        for(k=0; k < K; k++) {
            alphaHat_t(k, 0) = alphaStart(k);
        }
    }
    c(0) = 1 / sum(alphaHat_t.col(0));
    alphaHat_t.col(0) = alphaHat_t.col(0) * c(0);
    //
    // t here is 0-based
    //
    for(int t = 1; t < T; t++) {
        alphaConst = transMatRate_t_H(1, t-1) * arma::sum(alphaHat_t.col(t - 1));
        //
        alphaHat_t.col(t) = eMatHapSNP_t.col(t) % ( \
            transMatRate_t_H(0, t - 1) * alphaHat_t.col(t - 1) + \
            alphaConst * alphaMat_t.col(t - 1) );
        //
        c(t) = 1 / arma::sum(alphaHat_t.col(t));
        alphaHat_t.col(t) *= c(t);
    }
    return ;
}


//' @export
// [[Rcpp::export]]
void Rcpp_run_backward_haploid(
    arma::mat& betaHat_t,
    arma::rowvec& c,
    const arma::mat& eMatHapSNP_t,
    const arma::mat& alphaMat_t,    
    const arma::mat& transMatRate_t_H
) {
    const int T = eMatHapSNP_t.n_cols;
    //const int K = eMatHapSNP_t.n_rows;
    double x;
    arma::colvec e_times_b;
    for(int t = T-2; t >= 0; --t) {
        //x = 0;
        e_times_b = eMatHapSNP_t.col(t+1) % betaHat_t.col(t+1);
        x = transMatRate_t_H(1, t) * sum(alphaMat_t.col(t) % e_times_b);
        betaHat_t.col(t) = c(t) * (x + transMatRate_t_H(0, t) * e_times_b);
    }
    return;
}




//' @export
// [[Rcpp::export]]
arma::mat rcpp_make_eMatHap_t(
    const Rcpp::List& sampleReads,
    const int nReads,
    const arma::mat& eHaps_t,
    const double maxDifferenceBetweenReads,
    const int Jmax,
    arma::mat& eMatHapOri_t,
    const arma::vec& pRgivenH1,
    const arma::vec& pRgivenH2,
    const bool run_pseudo_haploid = false
) {
    //
    // constants
    //
    //const int T = eHaps_t.n_cols;
    const int K = eHaps_t.n_rows; // traditional K for haplotypes
    //
    // new variables
    //
    double pR = 0;
    double pA = 0;
    double eps, x;
    int j, k, J, readSNP, jj;
    arma::mat eMatHap_t = arma::ones(K,nReads);
    //
    // now build
    //
    int iRead;
    for(iRead=0; iRead<=nReads-1; iRead++) {
        Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
        // recal that below is what is used to set each element of sampleRead
        // note - this is no longer quite accurate
        // sampleReads.push_back(Rcpp::List::create(nU,d,phiU,pRU));
        J = as<int>(readData[0]); // number of Unique SNPs on read
        readSNP = as<int>(readData[1]); // leading SNP from read
        arma::ivec bqU = as<arma::ivec>(readData[2]); // bq for each SNP
        arma::ivec pRU = as<arma::ivec>(readData[3]); // position of each SNP from 0 to T-1
        // once each SNP is done, have P(read | k), can multiply to get P(read|(k1,k2))
        if(J>=Jmax)
            J=Jmax;
        for(j=0; j<=J; j++) {
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
            jj=pRU(j);
            for(k=0; k<=K-1; k++) {
                eMatHap_t(k,iRead) = eMatHap_t(k,iRead) * \
                    ( eHaps_t(k, jj) * pA + (1-eHaps_t(k,jj)) * pR);
            }
            if (run_pseudo_haploid == true) {
                //
                // pseudo-haploid          
                //
                x=pRgivenH1(iRead) / (pRgivenH1(iRead) + pRgivenH2(iRead));
                //
                for(k=0; k<=K-1; k++) {
                    // inefficient?
                    eMatHapOri_t(k, iRead) = eMatHap_t(k, iRead);
                    eMatHap_t(k,iRead) = x * eMatHap_t(k,iRead) + (1-x) * pRgivenH2(iRead);
                }
            }
        } 
        //
        // cap P(read|k) to be within maxDifferenceBetweenReads orders of magnitude
        //
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
    return(eMatHap_t);
}
    


//' @export
// [[Rcpp::export]]
arma::mat rcpp_make_eMatHapSNP_t(
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
    //const Rcpp::IntegerVector& wif
    int nReads = sampleReads.size(); //
    const int K = eMatHap_t.n_rows; // traditional K for haplotypes        
    arma::mat eMatHapSNP_t = arma::ones(K, nGrids); // why is this called SNP? 
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
    return(eMatHapSNP_t);
}


//' @export
// [[Rcpp::export]]
void make_haploid_gammaUpdate_t(
    arma::cube& gammaUpdate_t,
    const Rcpp::List& sampleReads,
    const int nReads,
    const arma::mat& gamma_t,
    const arma::mat& eHapsCurrent_t,
    const arma::mat& eMatHap_t,    
    const arma::mat& eMatHapOri_t,
    const arma::vec& pRgivenH1,
    const arma::vec& pRgivenH2,
    const bool run_pseudo_haploid = false
) {
    //
    //const int nSNPs = eHapsCurrent_t.n_cols;
    const int K = eHapsCurrent_t.n_rows;
    //
    //arma::cube gammaUpdate_t = arma::zeros(K, nSNPs, 2);
    int iRead, J, cr, t, j, k;
    Rcpp::List readData;
    arma::ivec bqU, pRU;
    arma::colvec gamma_t_col;
    double d3, eps, a1, a2, y, d, d1, d2, b;
    double pR = -1;
    double pA = -1;
    int cr_prev = -1;
    //
    for(iRead = 0; iRead < nReads; iRead++) {
        readData = as<Rcpp::List>(sampleReads[iRead]);
        J = as<int>(readData[0]); // number of SNPs on read
        cr = as<int>(readData[1]); // central SNP in read
        bqU = as<arma::ivec>(readData[2]); // bq for each SNP
        pRU = as<arma::ivec>(readData[3]); // position of each SNP from 0 to T-1
        if (run_pseudo_haploid == true) {      
          d3 = pRgivenH1(iRead) / (pRgivenH1(iRead) + pRgivenH2(iRead));
        }
        if (cr > cr_prev) {
            // do not always need to update
            gamma_t_col = gamma_t.col(cr);
        }
        for(j=0; j<=J; j++) {
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
            if (run_pseudo_haploid == true) {
                for(k=0; k<=K-1;k++) {
                    a1 = pA * eHapsCurrent_t(k,t);
                    a2 = pR * (1-eHapsCurrent_t(k,t)); 
                    y = a1 + a2;
                    b = eMatHapOri_t(k,iRead) * d3 / y;
                    d1 = a1 * b;
                    d2 = a2 * b;
                    d = gamma_t_col(k) / eMatHap_t(k,iRead);
                    gammaUpdate_t(k, t, 0) += d * d1;
                    gammaUpdate_t(k, t, 1) += d * (d1 + d2);
                }
            } else {
                for(k = 0; k < K; k++) {
                    a1 = pA * eHapsCurrent_t(k, t);
                    a2 = pR * (1 - eHapsCurrent_t(k, t));
                    y = a1 + a2;
                    d = gamma_t_col(k) / y;
                    gammaUpdate_t(k, t, 0) += a1 * d;
                    gammaUpdate_t(k, t, 1) += y * d;
                }
            }
        } // end of SNP in read 
    } // end of read
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
    









//' @export
// [[Rcpp::export]]
Rcpp::List forwardBackwardHaploid(
    const Rcpp::List& sampleReads,
    const int nReads,
    const arma::vec pi,
    const arma::mat& transMatRate_t_H,
    const arma::mat& alphaMat_t,
    const arma::mat& eHaps_t,
    arma::mat& alphaHat_t,
    arma::mat& betaHat_t,
    const double maxDifferenceBetweenReads,
    const double maxEmissionMatrixDifference,
    const int Jmax,
    const int suppressOutput,
    const int model,
    arma::cube& gammaUpdate_t, // see use_supplied_gammaUpdate_t
    arma::mat& jUpdate_t,
    arma::mat& hapSum_t,
    Rcpp::NumericVector& priorSum,
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
  const int T = alphaMat_t.n_cols + 1;  // what we iterate over / grid
  const int nGrids = T; // will make change eventually
  const int K = eHaps_t.n_rows; // traditional K for haplotypes
  //
  // new variables
  //
  // variables working on nReads
  if (!pass_in_alphaBeta) {
      alphaHat_t = arma::zeros(K, T);
      betaHat_t = arma::zeros(K, T);
  }
  arma::rowvec c = arma::zeros(1,T);  
  // eMatHapSNP works on the SNPs themselves
  arma::mat gamma_t = arma::zeros(K, T);  
  // variables for transition matrix and initialization
  // int variables and such
  int t, k;
  Rcpp::List alphaBetaBlocks;
  Rcpp::List to_return;  
  //
  //
  //
  // eMatHap
  //
  //
  //
  next_section="Make eMatHap";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  //
  // define eMatHap to work on the number of reads
  arma::mat eMatHapOri_t;
  if (run_pseudo_haploid) {
      eMatHapOri_t = arma::zeros(K, nReads);
  }
  arma::mat eMatHap_t = rcpp_make_eMatHap_t(
      sampleReads,
      nReads,
      eHaps_t,
      maxDifferenceBetweenReads,
      Jmax,
      eMatHapOri_t,
      pRgivenH1,
      pRgivenH2,
      run_pseudo_haploid
  );
  //
  //
  // once we have eMatHap, ie probabilities from reads, make eMatHapSNPsfrom this
  //
  //
  next_section="Make eMat";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  arma::mat eMatHapSNP_t = rcpp_make_eMatHapSNP_t(
      eMatHap_t,
      1,
      sampleReads,
      1,
      nGrids,
      run_fb_grid_offset,
      true,
      true,
      maxEmissionMatrixDifference,
      rescale
  );
  //
  //
  // forward recursion
  //
  //
  next_section="Forward recursion";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  //
  Rcpp_run_forward_haploid(
      alphaHat_t,
      c,
      eMatHapSNP_t,
      alphaMat_t,
      transMatRate_t_H,
      T,
      K,
      pi,
      alphaStart,
      run_fb_subset
  );
  //
  //
  // backward recursion
  //
  //
  next_section="Backward recursion";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  if (run_fb_subset == false) {    
      betaHat_t.col(T-1).fill(c(T-1));
  } else {
      for(k=0; k<=K-1; k++) {
          betaHat_t(k, T-1) = betaEnd(k);
      }
  }
  Rcpp_run_backward_haploid(
      betaHat_t,
      c,
      eMatHapSNP_t,
      alphaMat_t,
      transMatRate_t_H
  );
  //
  //
  // make gamma
  //
  //
  next_section="Make gamma";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  //
  gamma_t = alphaHat_t % betaHat_t;
  // normalize as well
  double g_temp = 0;
  for(t = 0; t < T; t++) {
      g_temp = 1 / c(t);
      gamma_t.col(t) *= g_temp;
  }
  //
  // optional early return, for final iteration
  //
  if (output_haplotype_dosages) {
      // calculate expanded version? for output?
      arma::mat gammaEK_t = make_gammaEK_t_from_gammaK_t(
          gamma_t, K, grid,
          snp_start_1_based, snp_end_1_based, run_fb_grid_offset
      );
      to_return.push_back(gammaEK_t, "gammaEK_t");      
  }
  //
  to_return.push_back(gamma_t, "gamma_t");
  //
  if (run_fb_subset == true) {
      return(to_return);
  }
  // make outputs here
  if (generate_fb_snp_offsets == true) {
      alphaBetaBlocks = rcpp_make_fb_snp_offsets(
          alphaHat_t,
          betaHat_t,
          blocks_for_output
      );
  }
  if (update_in_place) {  
      hapSum_t += gamma_t; // same as gammaK_t
      for(k=0; k < K; k++) {
          priorSum(k) = priorSum(k) + gamma_t(k, 0);
      }
  }
  //
  //
  //
  // hap probs are gamma!
  //
  //
  next_section="Gamma update";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  //
  if (!update_in_place) {        
      const int nSNPs = eHaps_t.n_cols;
      gammaUpdate_t = arma::zeros(K, nSNPs, 2);
  }
  make_haploid_gammaUpdate_t(  
      gammaUpdate_t,
      sampleReads,
      nReads,
      gamma_t,
      eHaps_t,
      eMatHap_t,      
      eMatHapOri_t,
      pRgivenH1,
      pRgivenH2,
      run_pseudo_haploid
  );
  //
  //
  // make jUpdate
  //
  //
  next_section="Make jUpdate";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  //
  if (!update_in_place) {
      jUpdate_t = arma::zeros(K, nGrids - 1);
  }
  for(t = 0; t < T - 1; t++) {
      jUpdate_t.col(t) += transMatRate_t_H(1, t) * (alphaMat_t.col(t) % betaHat_t.col(t + 1) % eMatHapSNP_t.col(t + 1));
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
      to_return.push_back(jUpdate_t, "jUpdate_t");
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
Rcpp::List rcpp_sample_multiple_paths(const int n_starts, const int n_its, const Rcpp::List& sampleReads, const int nReads, const arma::mat& eHaps_t, const double maxDifferenceBetweenReads, const int Jmax, const arma::vec pi, const arma::mat& transMatRate_t, const arma::mat& alphaMat_t, const arma::ivec& srp, const arma::ivec& sum_dosage_vec) {
    //
    // new variables go here
    //
    double n_save_iterations;
    const int T = eHaps_t.n_cols;
    int i_start, iRead, iHap, it, t;
    double s1, s2, p1, p2, pHap1;    
    arma::imat path(2, T);
    arma::mat read_labels = arma::zeros(2, nReads);
    arma::mat dosages = arma::zeros(n_starts, T);
    //
    // output variables go here
    //
    arma::mat eMatHapPH_t;
    arma::vec pRgivenH1;
    arma::vec pRgivenH2;
    arma::mat eMatHap_t = rcpp_make_eMatHap_t(
        sampleReads,
        nReads,
        eHaps_t,
        maxDifferenceBetweenReads,
        Jmax,
        eMatHapPH_t,
        pRgivenH1,
        pRgivenH2
    );
    //
    // initialize random numbers
    //
    double rand_uniform = 0;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    //
    // loop over starts
    //
    for(i_start = 0; i_start < n_starts; i_start++) {
        // start with random labels
        for(iRead=0; iRead < nReads; iRead++) {
            rand_uniform = dis(gen);
            if (rand_uniform > 0.5) {
                read_labels(0, iRead) = 1;
            } else {
                read_labels(1, iRead) = 1;
            }
        }
        for(it = 0; it < n_its; it++) {
            // sample path
            for(iHap = 0; iHap < 2; iHap++) {
                path.row(iHap) = rcpp_sample_path(read_labels.row(iHap), eMatHap_t, sampleReads, nReads, eHaps_t, maxDifferenceBetweenReads, Jmax, pi, transMatRate_t, alphaMat_t);
            }
            // sample labels here
            for(iRead = 0; iRead < nReads; iRead ++) {
                s1 = path(0, srp(iRead));
                s2 = path(1, srp(iRead));
                p1 = 0.5 * eMatHap_t(s1, iRead);
                p2 = 0.5 * eMatHap_t(s2, iRead);
                pHap1 = p1 / (p1 + p2);
                rand_uniform = dis(gen);
                if (rand_uniform < pHap1) {
                    read_labels(0, iRead) = 1;
                    read_labels(1, iRead) = 0;
                } else {
                    read_labels(0, iRead) = 0;
                    read_labels(1, iRead) = 1;
                }
            }
            // can save dosages here
            if (sum_dosage_vec(it) == 1) {
                for(iHap = 0; iHap < 2; iHap++) {
                    for(t=0; t < T; t++) {
                        dosages(i_start, t) = dosages(i_start, t) +      \
                            eHaps_t(path(iHap, t), t);
                    }
                }
            }

        } // end of iterations
        //
        // normalize dosage here
        //
        n_save_iterations = 0;
        for(int i=0; i < n_its; i++)
            n_save_iterations = n_save_iterations + sum_dosage_vec(i);
        for(t=0; t < T; t++) {
            dosages(i_start, t) = dosages(i_start, t) / n_save_iterations;
        }
    } // end of starts
    return(wrap(Rcpp::List::create(
                                   Rcpp::Named("dosages") = dosages,
                                   Rcpp::Named("path") = path,
                                   Rcpp::Named("read_labels") = read_labels
                                   )));
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

