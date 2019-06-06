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
    double& prev,
    const int suppressOutput,
    std::string& prev_section,
    std::string& next_section,
    const int grid_offset = 0
);



//' @export
// [[Rcpp::export]]
void collapse_diploid_gamma(
    arma::mat& gamma_t,
    arma::mat& gammaK_t,
    double& prev,
    const int suppressOutput,
    std::string& prev_section,
    std::string& next_section
) {
    next_section="collapse diploid gamma";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    int K_times_k1;
    int t, k1, k2;
    double d;
    const int nGrids = gamma_t.n_cols;
    const int K = gammaK_t.n_rows;
    for(t = 0; t < nGrids; t++) {
        for(k1 = 0; k1 < K; k1++) {
            d = 0;
            K_times_k1 = K * k1;
            for(k2 = 0; k2 < K; k2++) {
                d += gamma_t(K_times_k1 + k2, t);
            }
            gammaK_t(k1, t) = d;
        }
    }
    return;
}


//' @export
// [[Rcpp::export]]
void rcpp_make_and_bound_eMatGrid_diploid_t(
    arma::mat& eMatGrid_t,
    const arma::mat& eMatRead_t,
    const Rcpp::List& sampleReads,
    const double& maxEmissionMatrixDifference,
    double& prev,
    const int suppressOutput,
    std::string& prev_section,
    std::string& next_section,
    const int run_fb_grid_offset = 0,
    const bool rescale_eMatGrid_t = true
) {
    next_section="Initialize and bound eMatGrid";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //
    const int nReads = sampleReads.size(); //
    const int K = eMatRead_t.n_rows; // traditional K for haplotypes
    const int KK = K * K;
    const int nGrids = eMatGrid_t.n_cols;
    //
    int readSNP;
    double x, rescale;
    arma::colvec eMatRead_t_col;
    int iGrid, k3, k, k1, k2, iRead;
    for(iRead=0; iRead < nReads; iRead++) {
        Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
        readSNP = as<int>(readData[1]) - run_fb_grid_offset; // leading SNP from read
        eMatRead_t_col = 0.5 * eMatRead_t.col(iRead);
        for(k1 = 0; k1 < K; k1++) {
            x = eMatRead_t_col(k1);
            k3 = K * k1;
            for(k2 = 0; k2 < K; k2++) {
                eMatGrid_t(k3 + k2, readSNP) *= (x + eMatRead_t_col(k2));                
            }
        }
    }
    //
    if (!rescale_eMatGrid_t) {
        return;
    }
    //
    // cap eMat, ie P(reads | k1,k2) to be within maxDifferenceBetweenReads^2
    //
    double one_over_maxEmissionMatrixDifference = 1 / maxEmissionMatrixDifference;
    // loop over eMatGrid_t
    for(iGrid = 0; iGrid < nGrids; iGrid++) {
        if (eMatGrid_t(0, iGrid) < 1) {
            // first, get maximum
            x = eMatGrid_t.col(iGrid).max();
            // x is the maximum now. re-scale to x
            rescale = 1 / x;        
            for(k=0; k < KK; k++) {
                eMatGrid_t(k, iGrid) *= rescale;
                if(eMatGrid_t(k, iGrid)<(one_over_maxEmissionMatrixDifference)) {
                    eMatGrid_t(k, iGrid) = one_over_maxEmissionMatrixDifference;
                }
            }
        }
    } // end of loop on t
    return;
}








//' @export
// [[Rcpp::export]]
arma::imat sample_diploid_path(const arma::mat & alphaHat_t, const arma::mat & transMatRate_t_D, const arma::mat & eMatGrid_t, const arma::mat & alphaMat_t, const int T, const int K, const arma::rowvec & c) {
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
        norm = alphaHat_t(prev_state, t + 1) / eMatGrid_t(prev_state, t + 1) / c(t + 1);
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
    arma::cube& alphaMatSum_tc,
    int s,
    const arma::mat& alphaHat_t,
    const arma::mat& betaHat_t,    
    const arma::cube& transMatRate_tc_D,
    const arma::cube& alphaMatCurrent_tc,
    const arma::mat& eMatGrid_t,
    double& prev,
    const int suppressOutput,
    std::string& prev_section,
    std::string& next_section
) {
    next_section="Make xi-like calculations";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //
    const int nGrids = eMatGrid_t.n_cols;
    const int K = alphaMatCurrent_tc.n_rows; // traditional K for haplotypes
    //
    int iGrid, k1, k2, kk, K_times_k1;
    arma::vec alphaTemp1 = arma::zeros(K);
    arma::vec alphaTemp2 = arma::zeros(K);
    double tmr1, tmr2, d;
    arma::colvec alphaMat_t_col, alphaMat_t_col_times_tmr2, betaHat_times_eMat;
    //
    for(iGrid = 0; iGrid < nGrids - 1; iGrid++) {
        alphaTemp1.fill(0);
        for(k1 = 0; k1 < K; k1++) {
            K_times_k1 = K * k1;
            for(k2 = 0; k2 < K; k2++) {
                alphaTemp1(k2) += alphaHat_t(K_times_k1 + k2, iGrid);
            }
        }
        //
        // now do proper calculation
        //
        tmr1 = transMatRate_tc_D(1, iGrid, s);
        tmr2 = transMatRate_tc_D(2, iGrid, s);
        alphaMat_t_col = alphaMatCurrent_tc.slice(s).col(iGrid);
        alphaMat_t_col_times_tmr2 = alphaMat_t_col * tmr2;
        alphaTemp1 *= tmr1;
        alphaTemp1 += alphaMat_t_col_times_tmr2;
        betaHat_times_eMat = betaHat_t.col(iGrid + 1) % eMatGrid_t.col(iGrid + 1);
        for(k1 = 0; k1 < K; k1++) {
            K_times_k1 = K * k1;
            d = 0;
            for(k2 = 0; k2 < K; k2++) {
                kk = K_times_k1 + k2;
                d += alphaTemp1(k2) * betaHat_times_eMat(kk);
                // jUpdate_t(k1, t) += alphaTemp1(k2) * betaHat_t(kk, t + 1) * eMatGrid_t(kk, t + 1);
            }
            d *= 2 * alphaMat_t_col(k1);
            alphaMatSum_tc(k1, iGrid, s) += d;
        }   // end of loop on k
    } // end of loop on t
    return;
};




// requires initialization of first column
void run_forward_diploid(
    arma::mat& alphaHat_t,
    const Rcpp::NumericVector alphaStart,    
    arma::rowvec& c,
    const arma::mat& eMatGrid_t,
    const arma::cube& alphaMatCurrent_tc,
    const arma::cube& transMatRate_tc_D,
    const arma::mat& priorCurrent_m,
    int s,
    bool run_fb_subset,
    double& prev,
    int suppressOutput,
    std::string& prev_section,
    std::string& next_section
) {
    next_section="Forward recursion";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //
    const int nGrids = alphaMatCurrent_tc.n_cols + 1;
    const int K = alphaMatCurrent_tc.n_rows; // traditional K for haplotypes
    const int KK = K*K; // KK is number of states / traditional K for HMMs
    int kk, k1, k2, K_times_k1, iGrid;
    double alphaConst;
    arma::vec alphaTemp1 = arma::zeros(K);
    arma::vec alphaTemp2 = arma::zeros(K);
    arma::colvec alphaHat_t_col, alphaMat_t_col;
    double d0, d1, d2;
    //
    if (!run_fb_subset) {
        for(k1 = 0; k1 < K; k1++) {
            for(k2 = 0; k2 < K; k2++) {
                alphaHat_t(k1 + K * k2,0) = priorCurrent_m(k1, s) * priorCurrent_m(k2, s) * eMatGrid_t(k1 + K * k2, 0);
            }
        }
    } else {
        for(kk = 0; kk < KK; kk++) {
            alphaHat_t(kk, 0) = alphaStart(kk);
        }
    }
    c(0) = 1 / sum(alphaHat_t.col(0));
    alphaHat_t.col(0) = alphaHat_t.col(0) * c(0);
    //
    for(iGrid = 1; iGrid < nGrids; iGrid++) {
        // calculate necessary things
        alphaHat_t_col = alphaHat_t.col(iGrid - 1);
        alphaMat_t_col = alphaMatCurrent_tc.slice(s).col(iGrid - 1);
        d0 = transMatRate_tc_D(0, iGrid - 1, s);
        d1 = transMatRate_tc_D(1, iGrid - 1, s);
        d2 = transMatRate_tc_D(2, iGrid - 1, s);
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
                alphaHat_t(kk, iGrid) = eMatGrid_t(kk, iGrid) * \
                    (alphaHat_t_col(kk) * d0 +                  \
                     (alphaMat_t_col(k1) * alphaTemp1(k2) +     \
                      alphaMat_t_col(k2) * alphaTemp2(k1)) +            \
                     alphaMat_t_col(k1) * alphaMat_t_col(k2) * alphaConst);
            }
        }
        // do scaling now
        c(iGrid) = 1 / sum(alphaHat_t.col(iGrid));
        alphaHat_t.col(iGrid) *= c(iGrid);
    }
    return;
}


// requires initialization of first column
void run_backward_diploid(
    arma::mat& betaHat_t,
    arma::rowvec& c,
    const arma::mat& eMatGrid_t,
    const arma::cube& alphaMatCurrent_tc,
    const arma::cube& transMatRate_tc_D,
    int s,
    bool run_fb_subset,
    const Rcpp::NumericVector& betaEnd,
    double& prev,
    const int suppressOutput,
    std::string& prev_section,
    std::string& next_section
) {
    //
    next_section="Backward recursion";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //
    const int K = alphaMatCurrent_tc.n_rows; // traditional K for haplotypes
    const int KK = K*K; // KK is number of states / traditional K for HMMs
    const int nGrids = alphaMatCurrent_tc.n_cols + 1;  // what we iterate over / grid
    //
    int iGrid, k1, k2, kk, K_times_k1;
    arma::vec betaTemp1 = arma::zeros(K);
    arma::vec betaTemp2 = arma::zeros(K);
    double betaConst, d, x, d0;
    arma::colvec alphaMat_t_col, betaHat_mult_eMatGrid_t_col;
    //
    if (!run_fb_subset) {  
          betaHat_t.col(nGrids - 1).fill(c(nGrids - 1));
    } else {
        for(kk = 0; kk < KK; kk++) {
            betaHat_t(kk, nGrids - 1) = betaEnd(kk);
        }
    }
    //
    for(iGrid = nGrids - 2; iGrid>=0; --iGrid) {
        betaTemp1.fill(0);
        betaTemp2.fill(0);
        betaConst=0;
        alphaMat_t_col = alphaMatCurrent_tc.slice(s).col(iGrid);
        betaHat_mult_eMatGrid_t_col = betaHat_t.col(iGrid + 1) % eMatGrid_t.col(iGrid + 1);
        for(k1 = 0; k1 < K; k1++) {
            K_times_k1 = K * k1;
            for(k2 = 0; k2 < K; k2++) {
                kk = K_times_k1 + k2;
                d = betaHat_mult_eMatGrid_t_col(kk);
                betaTemp1(k1) += d * alphaMat_t_col(k2);
                betaTemp2(k2) += d * alphaMat_t_col(k1);
                betaConst += d * alphaMat_t_col(k1) * alphaMat_t_col(k2);
            }
        }
        // add transMatRate to constants
        d = transMatRate_tc_D(1, iGrid, s);
        betaTemp1 *= d;
        betaTemp2 *= d;
        betaConst *= transMatRate_tc_D(2, iGrid, s);
        d0 = transMatRate_tc_D(0, iGrid, s);
        // final calculation
        for(k1 = 0; k1 < K; k1++) {
            x = betaTemp1(k1) + betaConst;
            K_times_k1 = K * k1;            
            for(k2 = 0; k2 < K; k2++) {
                kk = K_times_k1 + k2;                
                betaHat_t(kk, iGrid) = betaHat_mult_eMatGrid_t_col(kk) * d0 + betaTemp2(k2) + x;
            }
        }
        // apply scaling
        betaHat_t.col(iGrid) *= c(iGrid);
    }
    return;
}




void make_diploid_gamma(
    arma::mat& gamma_t,
    arma::mat& alphaHat_t,
    arma::mat& betaHat_t,
    arma::rowvec& c,
    double& prev,
    const int suppressOutput,
    std::string& prev_section,
    std::string& next_section
) {
    next_section="Make gamma";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //
    gamma_t = alphaHat_t % betaHat_t;
    double g_temp = 0;
    for(int t = 0; t < gamma_t.n_cols; t++) {
        g_temp = 1 / c(t);
        gamma_t.col(t) *= g_temp;
    }
    return;
}

void calculate_diploid_gammaUpdate(
    arma::cube& gammaSum0_tc,
    arma::cube& gammaSum1_tc,
    int s,
    const Rcpp::List& sampleReads,
    const arma::mat& gamma_t,
    const arma::cube& eHapsCurrent_tc,
    const arma::mat& eMatRead_t,
    double& prev,
    const int suppressOutput,
    std::string& prev_section,
    std::string& next_section
) {
    //
    next_section="Gamma update";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    //
    const int K = eHapsCurrent_tc.n_rows;
    const int nReads = sampleReads.size();    
    //
    arma::colvec eMatRead_t_col;
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
            // RECAL  eMatRead(iRead,k) = eMatRead(iRead,k) * ( eHaps(pRU(j),k) * pA + (1-eHaps(pRU(j),k)) * pR);
            // RECAL  eMat(readSNP,k1+K*k2) = eMat(readSNP,k1+K*k2) * (0.5 * eMatRead(iRead,k1) + 0.5 * eMatRead(iRead,k2));
            //
            eMatRead_t_col = eMatRead_t.col(iRead);
            eHapsCurrent_t_col = eHapsCurrent_tc.slice(s).col(t);
            for(k1 = 0; k1 < K; k1++) {
                d1 = pA * eHapsCurrent_t_col(k1);
                d2 = pR * (1 - eHapsCurrent_t_col(k1));
                d3 = d1 / (d1 + d2); // this is all I need
                a = eMatRead_t_col(k1);
                K_times_k1 = K * k1;
                val1 = 0;
                val2 = 0;
                for(k2 = 0; k2 < K; k2++) {
                    b = eMatRead_t_col(k2);
                    kk = K_times_k1 + k2;                    
                    e = gamma_t_col(kk) * ( a / (a + b));
                    val1 += e * d3;
                    val2 += e;
                }
                gammaSum0_tc(k1, t, s) += 2 * val1; // this way I save the *2 multiplication until the end
                gammaSum1_tc(k1, t, s) += 2 * val2;
            } // end of loop onk
        } // end of loop on SNP within read
    }
    return;
}

//' @export
// [[Rcpp::export]]
void rcpp_calculate_fbd_dosage(
    arma::mat& genProbs_t,
    const arma::cube& eHapsCurrent_tc,
    int s,
    const arma::mat& gamma_t,
    const Rcpp::IntegerVector& grid,
    const int snp_start_1_based,
    const int snp_end_1_based,
    double& prev,
    const int suppressOutput,
    std::string& prev_section,
    std::string& next_section,
    const int grid_offset = 0
) {
    //
    next_section="Make genProbs_t";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
    // basically, copy and paste, either using grid or not
    const int nSNPs = snp_end_1_based - snp_start_1_based + 1;
    const int K = eHapsCurrent_tc.n_rows;
    // new
    int iSNP, t, k1, k2, kk, K_times_k1, cur_grid;
    double a, b, one_minus_b, g, one_minus_a, g0, g1, g2;
    arma::colvec eHapsCurrent_t_col, gamma_t_col, one_minus_eHapsCurrent_t_col;
    int prev_grid = -1;
    // i_t is index from 0 to nSNPs + 1, controls where things go
    // t is the index in the whole set of SNPs
    // tt is the index in the grid
    /// grid_offset refers to in what grids we are running
    for(iSNP = 0; iSNP < nSNPs; iSNP++) {
        g0 = 0;
        g1 = 0;
        g2 = 0;
        t = iSNP + snp_start_1_based - 1;
        cur_grid = grid(t) - grid_offset;
        if (cur_grid > prev_grid) {
            gamma_t_col = gamma_t.col(cur_grid);
            prev_grid = cur_grid;
        }
        eHapsCurrent_t_col = eHapsCurrent_tc.slice(s).col(t);
        one_minus_eHapsCurrent_t_col = 1 - eHapsCurrent_t_col;
        for(k1 = 0; k1 <K; k1++) {
            a = eHapsCurrent_t_col(k1);
            one_minus_a = one_minus_eHapsCurrent_t_col(k1);
            K_times_k1 = K * k1;
            for(k2 = 0; k2 < K; k2++) {
                kk = K_times_k1 + k2; // does not matter which way I do this
                g = gamma_t_col(kk);
                b = eHapsCurrent_t_col(k2);
                one_minus_b = one_minus_eHapsCurrent_t_col(k2);
                g0 += g * one_minus_a * one_minus_b;
                g1 += g * (a * one_minus_b + one_minus_a * b);
                g2 += g * a * b;
            }
        }
        genProbs_t(0, iSNP) += g0;
        genProbs_t(1, iSNP) += g1;
        genProbs_t(2, iSNP) += g2;        
        //for(j=1;j<=2;j++) {
        //    dosage(i_t) = dosage(i_t) + j * genProbs_t(j, i_t);
        //}
    }
    return;
}




//' @export
// [[Rcpp::export]]
Rcpp::List forwardBackwardDiploid(
    const Rcpp::List& sampleReads,
    const arma::cube& eHapsCurrent_tc,
    const arma::cube& alphaMatCurrent_tc,
    const arma::cube& transMatRate_tc_D,
    const arma::mat& priorCurrent_m,
    arma::mat& alphaHat_t,
    arma::mat& betaHat_t,
    arma::mat& gamma_t,
    arma::mat& eMatGrid_t,
    const double maxDifferenceBetweenReads,
    const double maxEmissionMatrixDifference,
    const int Jmax,
    const int suppressOutput,
    const arma::mat& blocks_for_output,
    arma::cube& gammaSum0_tc,
    arma::cube& gammaSum1_tc,    
    arma::cube& alphaMatSum_tc,
    arma::cube& hapSum_tc,
    arma::mat& priorSum_m,
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
    const bool output_haplotype_dosages = false, // whether to output state probabilities
    const bool rescale_eMatGrid_t = true // whether to rescale emat to minimize underflow problems, or to avoid, so log(c) = P(O | \lambda) has meaning    
) {
  double prev=clock();
  std::string prev_section="Null";
  std::string next_section="Initialize variables";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  //
  // constants
  //
  const int nGrids = alphaMatCurrent_tc.n_cols + 1;
  const int nSNPs = eHapsCurrent_tc.n_cols;
  const int K = eHapsCurrent_tc.n_rows;
  const int KK = K * K;
  const int S = eHapsCurrent_tc.n_slices;
  const int nReads = sampleReads.size();
  //
  // new variables
  //
  // variables working on nReads
  if (!pass_in_alphaBeta) {
      alphaHat_t = arma::zeros(KK, nGrids);
      betaHat_t = arma::zeros(KK, nGrids);
      gamma_t = arma::zeros(KK, nGrids);
      eMatGrid_t = arma::ones(KK, nGrids);
  }
  arma::rowvec c = arma::zeros(1, nGrids);
  // variables for faster forward backward calculation  double alphaConst, betaConst;
  arma::vec alphaTemp1 = arma::zeros(nGrids);
  arma::vec alphaTemp2 = arma::zeros(nGrids);
  arma::mat eMatRead_t = arma::ones(K, nReads);
  arma::mat genProbs_t = arma::zeros(3, nSNPs);  
  //
  // variables working on full space
  Rcpp::List alphaBetaBlocks;
  Rcpp::List to_return;
  Rcpp::List list_of_gamma_t;  
  if (!update_in_place) {
      gammaSum0_tc = arma::zeros(K, nSNPs, S);
      gammaSum1_tc = arma::zeros(K, nSNPs, S);      
      alphaMatSum_tc = arma::zeros(K, nGrids - 1, S);      
      hapSum_tc = arma::zeros(K, nGrids, S);
      priorSum_m = arma::zeros(K, S);
  }
  //
  //
  // dummy variables
  arma::mat gammaK_t;  
  if (!run_fb_subset | output_haplotype_dosages) {
      gammaK_t = arma::zeros(K, nGrids);
  }
  arma::mat gammaEK_t;
  arma::mat eMatHapPH_t;
  arma::vec pRgivenH1;
  arma::vec pRgivenH2;
  //
  //
  for(int s = 0; s < S; s++) {
      //
      if ((s > 0) | pass_in_alphaBeta) {
          alphaHat_t.fill(0);
          betaHat_t.fill(0);
          eMatRead_t.fill(1);
          eMatGrid_t.fill(1);
      }
      //
      rcpp_make_eMatRead_t(eMatRead_t, sampleReads, eHapsCurrent_tc, s, maxDifferenceBetweenReads, Jmax, eMatHapPH_t, pRgivenH1, pRgivenH2, prev, suppressOutput, prev_section, next_section);
      //
      rcpp_make_and_bound_eMatGrid_diploid_t(eMatGrid_t, eMatRead_t, sampleReads, maxEmissionMatrixDifference, prev, suppressOutput, prev_section, next_section, run_fb_grid_offset, rescale_eMatGrid_t);
      //
      run_forward_diploid(alphaHat_t, alphaStart, c, eMatGrid_t, alphaMatCurrent_tc, transMatRate_tc_D, priorCurrent_m, s, run_fb_subset, prev, suppressOutput, prev_section, next_section);
      //
      run_backward_diploid(betaHat_t, c, eMatGrid_t, alphaMatCurrent_tc, transMatRate_tc_D, s, run_fb_subset, betaEnd, prev, suppressOutput, prev_section, next_section);
      //
      make_diploid_gamma(gamma_t, alphaHat_t, betaHat_t, c, prev, suppressOutput, prev_section, next_section);
      if (return_gamma) {
          list_of_gamma_t.push_back(gamma_t, "gamma_t");      
      }
      //
      if (return_genProbs) {
          // seems convenient to wrap
          rcpp_calculate_fbd_dosage(genProbs_t, eHapsCurrent_tc, s, gamma_t, grid, snp_start_1_based, snp_end_1_based, prev, suppressOutput, prev_section, next_section, run_fb_grid_offset);
      }
      // make outputs here
      if (generate_fb_snp_offsets) {
          alphaBetaBlocks = rcpp_make_fb_snp_offsets(alphaHat_t, betaHat_t, blocks_for_output);
      }
      //
      // make collapsed gamma here (when is this wanted)
      // skip if run_fb_dosage and we don't want output haplotype dosages
      if (!run_fb_subset | output_haplotype_dosages) {
          // currently (probably?) only supported for S=1
          collapse_diploid_gamma(gamma_t, gammaK_t, prev, suppressOutput, prev_section, next_section);
          if (output_haplotype_dosages) {
              gammaEK_t = make_gammaEK_t_from_gammaK_t(gammaK_t, K, grid, snp_start_1_based, snp_end_1_based, prev, suppressOutput, prev_section, next_section, run_fb_grid_offset);
          }
      }
      //
      // optional, end early
      //
      if (!run_fb_subset) {
          //
          next_section="Update hapSum and prior";
          prev=print_times(prev, suppressOutput, prev_section, next_section);
          prev_section=next_section;
          hapSum_tc.slice(s) += gammaK_t;
          priorSum_m.col(s) += gammaK_t.col(0);
          //
          calculate_diploid_gammaUpdate(gammaSum0_tc, gammaSum1_tc, s, sampleReads, gamma_t, eHapsCurrent_tc, eMatRead_t, prev, suppressOutput, prev_section, next_section);
          //
          rcpp_make_diploid_jUpdate(alphaMatSum_tc, s, alphaHat_t, betaHat_t, transMatRate_tc_D, alphaMatCurrent_tc, eMatGrid_t, prev, suppressOutput, prev_section, next_section);
          //
      }
      //arma::imat sampled_path_diploid_t;
      //if (return_a_sampled_path) {
      //    sampled_path_diploid_t = sample_diploid_path(alphaHat_t, transMatRate_t_D, eMatGrid_t, alphaMat_t, T, K, c);
      //}
  }
  //
  //
  next_section="Done";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  if (return_genProbs) {
      if (S > 1) {
          // divide by S, here always just add
          genProbs_t *= 1 / double(S);
      }
      to_return.push_back(genProbs_t, "genProbs_t");
  }
  if (output_haplotype_dosages) {
      to_return.push_back(gammaEK_t, "gammaEK_t");
  }
  if (run_fb_subset) {
      // have already added genProbs hopefully?
      return(to_return);
  }
  if (!update_in_place) {
      to_return.push_back(gammaSum0_tc, "gammaSum0_tc");      
      to_return.push_back(gammaSum1_tc, "gammaSum1_tc");
      to_return.push_back(alphaMatCurrent_tc, "alphaMatCurrent_tc");
      to_return.push_back(hapSum_tc, "hapSum_tc");
      to_return.push_back(priorSum_m, "priorSum_m");      
  }
  if (generate_fb_snp_offsets) {
      to_return.push_back(alphaBetaBlocks, "alphaBetaBlocks");
  }
  // rest largely for debugging
  if (return_extra) {
      to_return.push_back(alphaHat_t, "alphaHat_t");
      to_return.push_back(betaHat_t, "betaHat_t");
      to_return.push_back(eMatRead_t, "eMatRead_t");      
      to_return.push_back(eMatGrid_t, "eMatGrid_t");
      to_return.push_back(c, "c");
  }
  if (return_gamma) {
      to_return.push_back(list_of_gamma_t, "list_of_gamma_t");
  }
  // deprecated?
  //if (return_a_sampled_path) {
  //to_return.push_back(sampled_path_diploid_t, "sampled_path_diploid_t");
  //}
  return(to_return);
}



