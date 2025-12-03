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


#define ERROR_INT -1



// suppose we have a vector prob of some length, summing to 1, and some single uniform
// return the integer val such that prob[val] < u < prob[val + 1]
int rcpp_simple_sample(
    Rcpp::NumericVector probs,
    double u
) {
  double v = 0;
  for(int i=0; i < probs.length(); i++) {
      v+= probs(i);
      if (u < v) {
          return(i);
      }
    }
  return ERROR_INT;
}


// [[Rcpp::export]]
NumericMatrix rcpp_phase_sample_paths_method_3(
    arma::cube& alphaMatCurrent_tc,
    arma::cube& transMatRate_tc_H,
    NumericMatrix alphaHat_t,
    NumericMatrix betaHat_t,
    NumericMatrix gamma_t,
    NumericMatrix eMatGrid_t,
    NumericVector c,
    int phasing_n_votes,
    arma::cube unifs_tc
) {
    int K = sqrt(alphaHat_t.nrow());
    double Kd = double(K);
    int KK = K * K;
    int nGrids = alphaHat_t.ncol();
    NumericMatrix paths(phasing_n_votes, nGrids);
    NumericVector one_through_KK(K * K);
    NumericVector amc(K);
    int k;
    for(k = 0; k < (KK); k++) {
      one_through_KK(k) = k + 1;
    }
    for(int i_phasing_vote = 0; i_phasing_vote < phasing_n_votes; i_phasing_vote++){
        paths(i_phasing_vote,0) = rcpp_simple_sample(gamma_t(_, 0), unifs_tc(i_phasing_vote, 0, 0)) + 1;
	for(int t2 = 1; t2 < nGrids; t2++) {
	    double U = unifs_tc(i_phasing_vote, t2, 0);
	    int t1 = t2 - 1; // 0-based
	    double subtrans1 = transMatRate_tc_H(0, t1, 0);
	    double subtrans2 = transMatRate_tc_H(1, t1, 0);
	    for(int k = 0; k < K; k++) {
	      amc(k) = alphaMatCurrent_tc(k, t1, 0); // assuming S = 1 (1-based)
	    }
	    int K1 = paths(i_phasing_vote, t1) - 1; // 0-based
	    int kt1_1 = K1 % K; // 0-based
	    int kt1_2 = std::ceil( (K1 + 1) / Kd) - 1; // 0-based
	    // if (i_phasing_vote == 3 && t2 == 1) {
	    //   for(k = 0; k < KK; k++) {
	    // 	int kt1_1 = k % K; // 0-based
	    // 	int kt1_2 = std::ceil( (k + 1) / Kd) - 1; // 0-based
	    // 	std::cout << "k = " << k << ", kt1_1 = " << kt1_1 << ", kt1_2 = " << kt1_2 << std::endl;
	    //   }
	    // }
	    double transprobstay = (subtrans2 * amc[kt1_1] + subtrans1)*(subtrans2 * amc[kt1_2] + subtrans1);
	    double pstay = transprobstay * eMatGrid_t(K1, t2) * betaHat_t(K1, t2) * c(t1) / betaHat_t(K1, t1);
	    //if (i_phasing_vote == 22) {
	    //  std::cout << "t2 = " << t2 << ",U = " << U << ", pstay = " << pstay << std::endl;
	    //}
	    if(U < pstay) {
	        paths(i_phasing_vote, t2) = K1 + 1;
	    } else{
	        NumericVector newprobs(K*K);
	        for(int kt2_1 = 0; kt2_1 < K; kt2_1++) {
		  for(int kt2_2 = 0; kt2_2 < K; kt2_2++) {
		    int K2 = kt2_1 + K * kt2_2;
		    double jump1 = subtrans2 * amc[kt2_1];
		    if (kt1_1 == kt2_1) {
		      jump1 += subtrans1;
		    }
		    double jump2 = subtrans2 * amc[kt2_2];
		    if (kt1_2 == kt2_2){
		      jump2 += subtrans1;
		    }
		    double trans_value = jump1 * jump2;
		    newprobs[K2] = trans_value * eMatGrid_t(K2, t2) * betaHat_t(K2, t2) * c(t1) / betaHat_t(K1, t1);
		  }
		}
		newprobs(K1) = 0;
		double s = Rcpp::sum(newprobs);
		for(k = 0; k < KK; k++) {
		  newprobs(k) /= s;
		}
                paths(i_phasing_vote, t2) = rcpp_simple_sample(newprobs, unifs_tc(i_phasing_vote, t2, 1)) + 1;
	    }
	}
    }
    return(paths);
}

