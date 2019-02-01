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
using namespace Rcpp;


double print_times(
    double prev,
    int suppressOutput,
    std::string past_text,
    std::string next_text
);

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
);

arma::mat rcpp_calculate_fbd_dosage(
    const arma::mat& eHapsCurrent_t,
    const arma::mat& gamma_t,
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

arma::mat make_gammaEK_t_from_gammaK_t(
    const arma::mat& gammaK_t,
    const int K,
    const Rcpp::IntegerVector& grid,
    const int snp_start_1_based,
    const int snp_end_1_based,
    const int grid_offset = 0
);



//' @export
// [[Rcpp::export]]
arma::mat collapse_diploid_gamma(
    const arma::mat& gamma_t,
    const int T,
    const int K
) {
    arma::mat gammaK_t = arma::zeros(K, T);
    int K_times_k1;
    int t, k1, k2;
    double d;
    for(t = 0; t < T; t++) {
        for(k1 = 0; k1 < K; k1++) {
            d = 0;
            K_times_k1 = K * k1;
            for(k2 = 0; k2 < K; k2++) {
                d += gamma_t(K_times_k1 + k2, t);
            }
            gammaK_t(k1, t) = d;
        }
    }
    return(gammaK_t);
}


//' @export
// [[Rcpp::export]]
arma::mat rcpp_make_and_bound_eMat_t(
    const arma::mat& eMatHap_t,
    const Rcpp::List& sampleReads,
    const int& nReads,
    const int& K,
    const int& T,
    const double& maxEmissionMatrixDifference,
    const int run_fb_grid_offset = 0
) {
    int readSNP;
    double x, rescale;
    const int KK = K * K;
    arma::mat eMat_t = arma::ones(KK, T);
    arma::colvec eMatHap_t_col;
    int iGrid, k3, k, k1, k2;
    for(int iRead=0; iRead < nReads; iRead++) {
        Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
        readSNP = as<int>(readData[1]) - run_fb_grid_offset; // leading SNP from read
        eMatHap_t_col = 0.5 * eMatHap_t.col(iRead);
        for(k1 = 0; k1 < K; k1++) {
            x = eMatHap_t_col(k1);
            k3 = K * k1;
            for(k2 = 0; k2 < K; k2++) {
                eMat_t(k3 + k2, readSNP) *= (x + eMatHap_t_col(k2));                
            }
        }  // end of SNP in read
    }
    //
    // cap eMat, ie P(reads | k1,k2) to be within maxDifferenceBetweenReads^2
    //
    double one_over_maxEmissionMatrixDifference = 1 / maxEmissionMatrixDifference;
    // loop over eMat_t
    for(iGrid = 0; iGrid < T; iGrid++) {
        if (eMat_t(0, iGrid) < 1) {
            // first, get maximum
            x = eMat_t.col(iGrid).max();
            // x is the maximum now. re-scale to x
            rescale = 1 / x;        
            for(k=0; k < KK; k++) {
                eMat_t(k, iGrid) *= rescale;
                if(eMat_t(k, iGrid)<(one_over_maxEmissionMatrixDifference))
                    eMat_t(k, iGrid) = one_over_maxEmissionMatrixDifference;
            }
        }
    } // end of loop on t
    return eMat_t;
}








//' @export
// [[Rcpp::export]]
arma::imat sample_diploid_path(const arma::mat & alphaHat_t, const arma::mat & transMatRate_t_D, const arma::mat & eMat_t, const arma::mat & alphaMat_t, const int T, const int K, const arma::rowvec & c) {
    //
    arma::imat sampled_path_diploid_t(3, T);
    int sampled_state;
    int t;
    int first_k;
    int second_k;
    int KK = K * K;
    double sum = 0;
    double rand_uniform = 0;
    double samp_vector_sum;
    double norm;
    arma::rowvec samp_vector = arma::zeros(1, KK);
    int prev_state;
    int prev_first_k, prev_second_k, k1, kk;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    // somehow this makes dis(gen) give me random 0-1 doubles I think
    //
    // choose initial state
    //
    double check = 0;
    for(kk = 0; kk < KK; kk++)
        check = check + alphaHat_t(kk, T - 1);
    if ((pow(check - 1, 2)) > 1e-4) {
        std::cout << "BAD input assumption" << std::endl;
        return(sampled_path_diploid_t);
    }
    //
    //
    kk = 0;  
    rand_uniform = dis(gen);
    sum = 0;
    while(kk < KK) {
        sum = sum + alphaHat_t(kk, T - 1); // does this sum to 1???
        if (sum > rand_uniform) {
            sampled_state = kk;
            first_k = (sampled_state) % K;
            second_k = (sampled_state - first_k) / K;
            sampled_path_diploid_t(0, T - 1) = first_k; // 0 based
            sampled_path_diploid_t(1, T - 1) = second_k; // 0 based
            sampled_path_diploid_t(2, T - 1) = sampled_state; // 0 based
            kk = KK;
        }
        kk = kk + 1;    
    }
    //
    // cycle
    //
    //for(t = T - 2; t >= 0; t--) {
    for(t = T - 2; t >= 0; t--) {
        prev_first_k = sampled_path_diploid_t(0, t + 1); // 0-based
        prev_second_k = sampled_path_diploid_t(1, t + 1); // 0 based
        prev_state = sampled_path_diploid_t(2, t + 1); // 0 based
        samp_vector.fill(transMatRate_t_D(2, t) * alphaMat_t(prev_first_k, t) * alphaMat_t(prev_second_k, t));
        for(k1=0; k1<=K-1; k1++) {
            // switch on first, keep second
            samp_vector(k1 + K * prev_second_k) = samp_vector(k1 + K * prev_second_k) + transMatRate_t_D(1, t) * alphaMat_t(prev_first_k, t);
            // keep first, switch on second
            samp_vector(prev_first_k + K * k1) = samp_vector(prev_first_k + K * k1) + transMatRate_t_D(1, t) * alphaMat_t(prev_second_k, t);
        }
        samp_vector(prev_state)=samp_vector(prev_state) + transMatRate_t_D(0, t);
        //
        samp_vector_sum = 0;
        for(kk=0; kk<=KK-1; kk++) {
            samp_vector(kk)=samp_vector(kk) * alphaHat_t(kk, t);
            samp_vector_sum=samp_vector_sum + samp_vector(kk);
        }
        norm = alphaHat_t(prev_state, t + 1) / eMat_t(prev_state, t + 1) / c(t + 1);
        if ((pow(samp_vector_sum / norm - 1, 2)) > 1e-4) {
             std::cout << "BAD COUNT on t=" << t << std::endl;
             std::cout << "samp_vector_sum=" << samp_vector_sum << std::endl;
             std::cout << "norm=" << norm << std::endl;
             return(sampled_path_diploid_t);
        }
        // sample state using crazy C++
        kk = 0;  
        rand_uniform = dis(gen);
        sum = 0;
        while(kk < KK) {
            sum = sum + samp_vector(kk) / norm;
            if (sum > rand_uniform) {
                sampled_state = kk;
                kk = KK;
            }
            kk = kk + 1;    
        }
        // convert back here
        first_k = (sampled_state) % K;
        second_k = (sampled_state - first_k) / K;
        sampled_path_diploid_t(0, t) = first_k;
        sampled_path_diploid_t(1, t) = second_k;
        sampled_path_diploid_t(2, t) = sampled_state;
    }
    return(sampled_path_diploid_t);
}



//' @export
// [[Rcpp::export]]
void rcpp_make_diploid_jUpdate(
    arma::mat& jUpdate_t,
    const int K,
    const int T,
    const arma::mat& alphaHat_t,
    const arma::mat& betaHat_t,    
    const arma::mat& transMatRate_t_D,
    const arma::mat& alphaMat_t,
    const arma::mat& eMat_t
) {
    int t, k1, k2, kk, K_times_k1;
    arma::vec alphaTemp1 = arma::zeros(K);
    arma::vec alphaTemp2 = arma::zeros(K);
    double tmr1, tmr2, d;
    arma::colvec alphaMat_t_col_times_tmr2, betaHat_times_eMat;
    //
    for(t=0; t<=T-2; t++) {
        alphaTemp1.fill(0);
        for(k1 = 0; k1 < K; k1++) {
            K_times_k1 = K * k1;
            for(k2 = 0; k2 < K; k2++) {
                alphaTemp1(k2) += alphaHat_t(K_times_k1 + k2, t);
            }
        }
        //
        // now do proper calculation
        //
        tmr1 = transMatRate_t_D(1, t);
        tmr2 = transMatRate_t_D(2, t);
        alphaMat_t_col_times_tmr2 = alphaMat_t.col(t) * tmr2;
        alphaTemp1 *= tmr1;
        alphaTemp1 += alphaMat_t_col_times_tmr2;
        betaHat_times_eMat = betaHat_t.col(t + 1) % eMat_t.col(t + 1);
        for(k1 = 0; k1 < K; k1++) {
            K_times_k1 = K * k1;
            d = 0;
            for(k2 = 0; k2 < K; k2++) {
                kk = K_times_k1 + k2;
                d += alphaTemp1(k2) * betaHat_times_eMat(kk);
                // jUpdate_t(k1, t) += alphaTemp1(k2) * betaHat_t(kk, t + 1) * eMat_t(kk, t + 1);
            }
            d *= 2 * alphaMat_t(k1, t);
            jUpdate_t(k1, t) += d;
        }   // end of loop on k
    } // end of loop on t
    return;
};




// requires initialization of first column
void run_forward_diploid(
    arma::mat& alphaHat_t,
    arma::rowvec& c,
    const arma::mat& eMat_t,
    const arma::mat& alphaMat_t,    
    const arma::mat& transMatRate_t_D,
    const int& T,
    const int& K
) {
    double alphaConst;
    int kk, k1, k2, K_times_k1;
    arma::vec alphaTemp1 = arma::zeros(K);
    arma::vec alphaTemp2 = arma::zeros(K);
    arma::colvec alphaHat_t_col, alphaMat_t_col;
    double d0, d1, d2;
    for(int t = 1; t < T; t++) {
        // calculate necessary things
        alphaHat_t_col = alphaHat_t.col(t - 1);
        alphaMat_t_col = alphaMat_t.col(t - 1);        
        d0 = transMatRate_t_D(0, t-1);
        d1 = transMatRate_t_D(1, t-1);
        d2 = transMatRate_t_D(2, t-1);
        alphaTemp1.fill(0);
        alphaTemp2.fill(0);
        for(k1 = 0; k1 < K; k1++) {
            K_times_k1 = K * k1;
            for(k2 = 0; k2 < K; k2++) {
                kk = K_times_k1 + k2;
                alphaTemp1(k2) += alphaHat_t_col(kk);
                alphaTemp2(k1) += alphaHat_t_col(kk);
            }
        }
        // now make constant over whole thing
        alphaConst = arma::sum(alphaHat_t_col) * d2;
        alphaTemp1 *= d1;
        alphaTemp2 *= d1;        
        //
        for(k1 = 0; k1 < K; k1++) {
            K_times_k1 = K * k1;            
            for(k2 = 0; k2 < K; k2++) {
                kk = K_times_k1 + k2;
                alphaHat_t(kk, t) = eMat_t(kk, t) *  \
                    (alphaHat_t_col(kk) * d0 +       \
                     (alphaMat_t_col(k1) * alphaTemp1(k2) + \
                      alphaMat_t_col(k2) * alphaTemp2(k1)) + \
                     alphaMat_t_col(k1) * alphaMat_t_col(k2) * alphaConst);
            }
        }
        // do scaling now
        c(t) = 1 / sum(alphaHat_t.col(t));
        alphaHat_t.col(t) *= c(t);
    }
    return;
}


// requires initialization of first column
void run_backward_diploid(
    arma::mat& betaHat_t,
    arma::rowvec& c,
    const arma::mat& eMat_t,
    const arma::mat& alphaMat_t,    
    const arma::mat& transMatRate_t_D,
    const int& T,
    const int& K
) {
    int t, k1, k2, kk, K_times_k1;
    arma::vec betaTemp1 = arma::zeros(K);
    arma::vec betaTemp2 = arma::zeros(K);
    double betaConst, d, x, d0;
    arma::colvec alphaMat_t_col, betaHat_mult_eMat_t_col;
    for(t = T-2; t>=0; --t) {
        betaTemp1.fill(0);
        betaTemp2.fill(0);
        betaConst=0;
        alphaMat_t_col = alphaMat_t.col(t);
        betaHat_mult_eMat_t_col = betaHat_t.col(t + 1) % eMat_t.col(t + 1); // element-wise multiplication
        for(k1 = 0; k1 < K; k1++) {
            K_times_k1 = K * k1;
            for(k2 = 0; k2 < K; k2++) {
                kk = K_times_k1 + k2;
                d = betaHat_mult_eMat_t_col(kk);
                betaTemp1(k1) += d * alphaMat_t_col(k2);
                betaTemp2(k2) += d * alphaMat_t_col(k1);
                betaConst += d * alphaMat_t_col(k1) * alphaMat_t_col(k2);
            }
        }
        // add transMatRate to constants
        d = transMatRate_t_D(1,t);
        betaTemp1 *= d;
        betaTemp2 *= d;
        betaConst *= transMatRate_t_D(2, t);
        d0 = transMatRate_t_D(0, t);
        // final calculation
        for(k1 = 0; k1 < K; k1++) {
            x = betaTemp1(k1) + betaConst;
            K_times_k1 = K * k1;            
            for(k2 = 0; k2 < K; k2++) {
                kk = K_times_k1 + k2;                
                betaHat_t(kk, t) = betaHat_mult_eMat_t_col(kk) * d0 + betaTemp2(k2) + x;
            }
        }
        // apply scaling
        betaHat_t.col(t) *= c(t);
    }
    return;
}




void calculate_diploid_gammaUpdate(
    arma::cube& gammaUpdate_t,
    const Rcpp::List& sampleReads,
    const int nReads,
    const arma::mat& gamma_t,
    const arma::mat& eHapsCurrent_t,
    const arma::mat& eMatHap_t 
) {
    //
    const int K = eHapsCurrent_t.n_rows;
    //
    arma::colvec eMatHap_t_col;
    arma::colvec eHapsCurrent_t_col, gamma_t_col;
    arma::ivec bqU, pRU;    
    int J, cr, j, iRead, k1, k2, t, K_times_k1, kk;
    int cr_prev = -1;
    double eps, pA, pR, d1, d2, d3, a, b, e, val1, val2;
    //
    for(iRead = 0; iRead < nReads; iRead++) {
        // recal that below is what is used to set each element of sampleRead
        // sampleReads.push_back(Rcpp::List::create(nU,d,phiU,pRU));
        Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
        J = readData[0]; // number of SNPs on read
        cr = readData[1]; // central SNP or grid point
        bqU = as<arma::ivec>(readData[2]); // bq for each SNP
        pRU = as<arma::ivec>(readData[3]); // position of each SNP from 0 to T-1
        // loop over every SNP in the read
        if (cr > cr_prev) {
            gamma_t_col = gamma_t.col(cr);
            cr_prev = cr;
        }
        for(j = 0; j <= J; j++) {
            t=pRU(j); // position of this SNP in full T sized matrix
            //
            // first haplotype (ie (k,k1))
            //
            // less than 0 - reference base more likely
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
            // RECAL  eMatHap(iRead,k) = eMatHap(iRead,k) * ( eHaps(pRU(j),k) * pA + (1-eHaps(pRU(j),k)) * pR);
            // RECAL  eMat(readSNP,k1+K*k2) = eMat(readSNP,k1+K*k2) * (0.5 * eMatHap(iRead,k1) + 0.5 * eMatHap(iRead,k2));
            //
            eMatHap_t_col = eMatHap_t.col(iRead);
            eHapsCurrent_t_col = eHapsCurrent_t.col(t);
            for(k1 = 0; k1 < K; k1++) {
                d1 = pA * eHapsCurrent_t_col(k1);
                d2 = pR * (1 - eHapsCurrent_t_col(k1));
                d3 = d1 / (d1 + d2); // this is all I need
                a = eMatHap_t_col(k1);
                K_times_k1 = K * k1;
                val1 = 0;
                val2 = 0;
                for(k2 = 0; k2 < K; k2++) {
                    b = eMatHap_t_col(k2);
                    kk = K_times_k1 + k2;                    
                    e = gamma_t_col(kk) * ( a / (a + b));
                    val1 += e * d3;
                    val2 += e;
                }
                gammaUpdate_t(k1, t, 0) += 2 * val1; // this way I save the *2 multiplication until the end
                gammaUpdate_t(k1, t, 1) += 2 * val2;
            } // end of loop onk
        } // end of loop on SNP within read
    }
    return;
}


//' @export
// [[Rcpp::export]]
Rcpp::List forwardBackwardDiploid(
    const Rcpp::List& sampleReads,
    const int nReads,
    const arma::vec& pi,
    const arma::mat& transMatRate_t_D,
    const arma::mat& alphaMat_t,
    const arma::mat& eHaps_t,
    arma::mat& alphaHat_t,
    arma::mat& betaHat_t,
    const double maxDifferenceBetweenReads,
    const double maxEmissionMatrixDifference,
    const int Jmax,
    const int suppressOutput,
    const arma::mat& blocks_for_output,
    arma::cube& gammaUpdate_t, // see update_in_place
    arma::mat& jUpdate_t,
    arma::mat& hapSum_t,
    Rcpp::NumericVector& priorSum,
    const bool generate_fb_snp_offsets = false,
    const Rcpp::NumericVector alphaStart = 0,
    const Rcpp::NumericVector betaEnd = 0,
    const bool return_a_sampled_path = false,
    const bool run_fb_subset = false,
    const int run_fb_grid_offset = 0, // this is 0-based
    const bool return_genProbs = false, // the below are needed if we want genProbs. dosage trivial
    const int snp_start_1_based = -1,
    const int snp_end_1_based = -1,
    const Rcpp::IntegerVector grid = 0, // end of things needed for genProbs and dosage
    const bool return_gamma = false, // full gamma, K * K rows
    const bool return_extra = false, // whether to return stuff useful for debugging
    const bool update_in_place = false, // update directly into output variables
    const bool pass_in_alphaBeta = false, // whether to pass in pre-made alphaHat, betaHat
    const bool output_haplotype_dosages = false // whether to output state probabilities
) {
  double prev=clock();
  std::string prev_section="Null";
  std::string next_section="Initialize variables";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  //
  // constants
  //
  const int T = alphaMat_t.n_cols + 1;  // what we iterate over / grid
  const int nSNPs = eHaps_t.n_cols; // traditional K for haplotypes  
  const int K = eHaps_t.n_rows; // traditional K for haplotypes
  const int KK = K*K; // KK is number of states / traditional K for HMMs
  //
  // new variables
  //
  // variables working on nReads
  if (!pass_in_alphaBeta) {
      alphaHat_t = arma::zeros(KK, T);  
      betaHat_t = arma::zeros(KK, T);
  }
  arma::rowvec c = arma::zeros(1, T);
  // variables for faster forward backward calculation  double alphaConst, betaConst;
  arma::vec alphaTemp1 = arma::zeros(K);
  arma::vec alphaTemp2 = arma::zeros(K);
  // variables working on full space
  int kk, t, k1, k2;
  Rcpp::List alphaBetaBlocks;
  Rcpp::List to_return;
  arma::mat genProbs_t;  
  //
  //
  // eMat - make complete
  //
  //
  next_section="Initialize eMatHap";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  // dummy variables
  arma::mat gammaK_t;
  arma::mat gammaEK_t;
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
  // once we have all the eMatHaps, ie probabilities from reads, make eMat from this
  //
  next_section="Initialize and bound eMat";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  //
  arma::mat eMat_t = rcpp_make_and_bound_eMat_t(eMatHap_t, sampleReads, nReads, K, T, maxEmissionMatrixDifference, run_fb_grid_offset);
  //
  // forward recursion
  //
  next_section="Forward recursion";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  if (run_fb_subset == false) {
      for(k1=0; k1<=K-1; k1++)
          for(k2=0; k2<=K-1; k2++)
              alphaHat_t(k1+K*k2,0) = pi(k1) * pi(k2) * eMat_t(k1+K*k2,0);
  } else {
      for(kk=0; kk<=KK-1; kk++) {
          alphaHat_t(kk, 0) = alphaStart(kk);
      }
  }
  c(0) = 1 / sum(alphaHat_t.col(0));
  alphaHat_t.col(0) = alphaHat_t.col(0) * c(0);
  run_forward_diploid(alphaHat_t, c, eMat_t, alphaMat_t, transMatRate_t_D, T, K);
  //
  //
  //
  // backward recursion
  //
  //
  next_section="Backward recursion";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  //
  if (run_fb_subset == false) {  
      betaHat_t.col(T-1).fill(c(T-1));
  } else {
      for(kk=0; kk<=KK-1; kk++) {
          betaHat_t(kk, T-1) = betaEnd(kk);
      }
  }
  run_backward_diploid(betaHat_t, c, eMat_t, alphaMat_t, transMatRate_t_D, T, K);
  //
  //
  // make gamma
  //
  //
  next_section="Make gamma";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  //
  arma::mat gamma_t = alphaHat_t % betaHat_t;
  double g_temp = 0;
  for(t = 0; t < T; t++) {
      g_temp = 1 / c(t);
      gamma_t.col(t) *= g_temp;
  }
  //
  // (optional) calculate genProbs and dosage
  //
  if (return_genProbs) {
      next_section="Make genProbs_t";
      prev=print_times(prev, suppressOutput, prev_section, next_section);
      prev_section=next_section;
      genProbs_t = rcpp_calculate_fbd_dosage(
          eHaps_t,
          gamma_t,
          grid,
          snp_start_1_based,
          snp_end_1_based,
          run_fb_grid_offset
      );
      to_return.push_back(genProbs_t, "genProbs_t");      
  }
  // make outputs here
  if (generate_fb_snp_offsets) {
      alphaBetaBlocks = rcpp_make_fb_snp_offsets(
          alphaHat_t,
          betaHat_t,
          blocks_for_output
      );
  }
  //
  // make collapsed gamma here
  //
  // skip if run_fb_dosage and we don't want output haplotype dosages
  if (!(run_fb_subset & !output_haplotype_dosages)) {
      next_section="Make collapsed gamma";
      prev=print_times(prev, suppressOutput, prev_section, next_section);
      prev_section=next_section;
      gammaK_t = collapse_diploid_gamma(gamma_t, T, K);
      if (output_haplotype_dosages) {
          gammaEK_t = make_gammaEK_t_from_gammaK_t(
              gammaK_t, K, grid,
              snp_start_1_based, snp_end_1_based, run_fb_grid_offset
          );
      }
  }
  //
  // optional, end early
  //
  if (run_fb_subset) {
      if (output_haplotype_dosages) {
          to_return.push_back(gammaEK_t, "gammaEK_t");
      }
      return(to_return);
  }
  //
  //
  //
  if (update_in_place) {  
      hapSum_t += gammaK_t;
      for(int k=0; k < K; k++) {
          priorSum(k) = priorSum(k) + gammaK_t(k, 0);
      }
  }
  //
  //
  // do gamma update here - save large matrix, just put in necessary values
  //
  //
  next_section="Gamma update";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  if (!update_in_place) {    
      gammaUpdate_t = arma::zeros(K, nSNPs, 2);      
  }
  calculate_diploid_gammaUpdate(
      gammaUpdate_t,
      sampleReads,
      nReads,
      gamma_t,
      eHaps_t,
      eMatHap_t 
  );
  //
  //
  // make jUpdate
  //
  //
  next_section="Make xi-like calculations";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  if (!update_in_place) {    
      jUpdate_t = arma::zeros(K, T - 1);
  }
  rcpp_make_diploid_jUpdate(
      jUpdate_t,
      K,
      T,
      alphaHat_t,
      betaHat_t,      
      transMatRate_t_D,
      alphaMat_t,
      eMat_t
  );
  //
  // optional, sample a path
  //
  arma::imat sampled_path_diploid_t;
  if (return_a_sampled_path) {
      sampled_path_diploid_t = sample_diploid_path(alphaHat_t, transMatRate_t_D, eMat_t, alphaMat_t, T, K, c);
  }
  //
  //
  // done - return necessary information
  //
  //
  next_section="Done";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  to_return.push_back(gammaK_t, "gammaK_t");
  if (!update_in_place) {    
      to_return.push_back(gammaUpdate_t, "gammaUpdate_t");
      to_return.push_back(jUpdate_t, "jUpdate_t");      
  }
  if (generate_fb_snp_offsets) {
      to_return.push_back(alphaBetaBlocks, "alphaBetaBlocks");
  }
  // rest largely for debugging
  if (return_extra) {
      to_return.push_back(alphaHat_t, "alphaHat_t");
      to_return.push_back(betaHat_t, "betaHat_t");
      to_return.push_back(eMat_t, "eMat_t");
      to_return.push_back(eMatHap_t, "eMatHap_t");
      to_return.push_back(c, "c");
  }
  if (return_gamma) { 
      to_return.push_back(gamma_t, "gamma_t");
  }
  // deprecated?
  if (return_a_sampled_path) {
      to_return.push_back(sampled_path_diploid_t, "sampled_path_diploid_t");
  }
  if (output_haplotype_dosages) {
      to_return.push_back(gammaEK_t, "gammaEK_t");
  }
  return(to_return);
}



