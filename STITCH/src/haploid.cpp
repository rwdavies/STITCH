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



//' @export
// [[Rcpp::export]]
arma::mat rcpp_make_eMatHap_t(const Rcpp::List& sampleReads, const int nReads, const arma::mat& eHaps_t, const double maxDifferenceBetweenReads, const int Jmax) {
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
        for(k=0; k<=K-1; k++)
            eMatHap_t(k,iRead) = eMatHap_t(k,iRead) *   \
                ( eHaps_t(k, jj) * pA + (1-eHaps_t(k,jj)) * pR);
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
arma::rowvec rcpp_sample_path(const arma::rowvec read_labels, const arma::mat eMatHap_t, const Rcpp::List& sampleReads, const int nReads, const arma::mat& eHaps_t, const double maxDifferenceBetweenReads, const int Jmax, const arma::vec pi, const arma::mat& transMatRate_t, const arma::mat& alphaMat_t) {
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
  arma::rowvec sampled_state = arma::zeros(1, T);
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
      alphaConst = alphaConst * transMatRate_t(1, t-1);
      //
      // each entry is emission * (no change * that value + constant)
      for(k=0; k<=K-1; k++)
          alphaHat_t(k,t) =  eMatHapSNP_t(k,t) *  \
              ( transMatRate_t(0, t-1) * alphaHat_t(k,t-1) +   \
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
      if ((pow(samp_vector_sum / norm, 2) - 1) > 1e-4) {
          std::cout << "BAD COUNT on t=" << t << std::endl;
          std::cout << "samp_vector_sum=" << samp_vector_sum << std::endl;
          std::cout << "norm=" << norm << std::endl;
          return(sampled_state);          
      }
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
Rcpp::List rcpp_sample_multiple_paths(const int n_starts, const int n_its, const Rcpp::List& sampleReads, const int nReads, const arma::mat& eHaps_t, const double maxDifferenceBetweenReads, const int Jmax, const arma::vec pi, const arma::mat& transMatRate_t, const arma::mat& alphaMat_t, const arma::ivec& srp, const arma::ivec& sum_dosage_vec) {
    //
    // new variables go here
    //
    double n_save_iterations;
    const int T = eHaps_t.n_cols;
    int i_start, iRead, iHap, it, t;
    double s1, s2, p1, p2, pHap1;    
    arma::mat path = arma::zeros(2, T);
    arma::mat read_labels = arma::zeros(2, nReads);
    arma::mat dosages = arma::zeros(n_starts, T);
    //
    // output variables go here
    //
    arma::mat eMatHap_t = rcpp_make_eMatHap_t(sampleReads, nReads, eHaps_t, maxDifferenceBetweenReads, Jmax);
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
