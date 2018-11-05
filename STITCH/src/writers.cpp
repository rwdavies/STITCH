#include <RcppArmadillo.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
Rcpp::StringVector rcpp_make_column_of_vcf(
    const arma::mat& gp_t,
    const bool use_read_proportions,
    const bool use_state_probabilities,
    const arma::mat& read_proportions,
    const arma::mat& q_t
) {
    // ## write out genotype, genotype likelihood, and dosage
    // ##GT:GL:DS
    // ##FORMAT=<ID=GT:,Number=1,Type=String,Description="Best Guessed Genotype with posterior probability threshold of 0.9">
    // ##FORMAT=<ID=GL,Number=3,Type=Float,Description="Posterior probability of 0/0, 0/1, and 1/1">
    // ##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">
    // ## 1/1:0,0.054,0.946:1.946
    // ## add one samples worth of info to a VCF
    const int T = gp_t.n_cols;
    Rcpp::StringVector output(T);
    int t = 0;
    int k;
    char buffer[30];
    const int K = q_t.n_rows;
    for(t=0; t < T; t++) {
      output(t) = "";
      // first, get genotype at the start
      if (gp_t(0, t) < gp_t(1, t)) {
	if (gp_t(1, t) < gp_t(2, t)) {
	  if (0.90 <= gp_t(2, t)) {
	    output(t) = "1/1";
	  } else {
	    output(t) = "./.";
	  }
	} else {
	  if (0.90 <= gp_t(1, t)) {
	    output(t) = "0/1";	  
	  } else {
	    output(t) = "./.";
	  }
	}
      } else {
	if (gp_t(0, t) < gp_t(2, t)) {
	  if (0.90 <= gp_t(2, t)) {
	    output(t) = "1/1";	  	  
	  } else {
	    output(t) = "./.";
	  }
	} else {
	  if (0.90 <= gp_t(0, t)) {
	    output(t) = "0/0";
	  } else {
	    output(t) = "./.";
	  }
	}
      }
      // now, add on gp
      sprintf(buffer, ":%.3f,%.3f,%.3f:%.3f", gp_t(0, t), gp_t(1, t), gp_t(2, t), gp_t(1, t) + 2 * gp_t(2, t));
      output(t) += buffer;
      if (use_read_proportions) {
	sprintf(buffer, ":%.3f,%.3f,%.3f,%.3f",
		read_proportions(t, 0),
		read_proportions(t, 1),
		read_proportions(t, 2),
		read_proportions(t, 3));
	output(t) += buffer;	
      }
      if (use_state_probabilities) {
	// ok, here, can I do each, one at a time?
	output(t) += ":";
	for(k = 0; k < K; k++) {
	  // add to string
	  sprintf(buffer, "%.3f", q_t(k, t));
	  output(t) += buffer;
	  // add comma is not last one
	  if (k < (K - 1)) {
	    output(t) += ",";	    
	  }
	}
      }
    }
    return output;
}
