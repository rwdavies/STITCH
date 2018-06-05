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
    const arma::mat& gp,
    const int use_read_proportions,
    const arma::mat& read_proportions
) {
    // ## write out genotype, genotype likelihood, and dosage
    // ##GT:GL:DS
    // ##FORMAT=<ID=GT:,Number=1,Type=String,Description="Best Guessed Genotype with posterior probability threshold of 0.9">
    // ##FORMAT=<ID=GL,Number=3,Type=Float,Description="Posterior probability of 0/0, 0/1, and 1/1">
    // ##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">
    // ## 1/1:0,0.054,0.946:1.946
    // ## add one samples worth of info to a VCF
    const int T = gp.n_rows;
    Rcpp::StringVector output(T);
    int t = 0;
    char buffer[30];
    for(t=0; t < T; t++) {
      output(t) = "";
      // first, get genotype at the start
      if (gp(t, 0) < gp(t, 1)) {
	if (gp(t, 1) < gp(t, 2)) {
	  if (0.90 <= gp(t, 2)) {
	    output(t) = "1/1";
	  } else {
	    output(t) = "./.";
	  }
	} else {
	  if (0.90 <= gp(t, 1)) {
	    output(t) = "0/1";	  
	  } else {
	    output(t) = "./.";
	  }
	}
      } else {
	if (gp(t, 0) < gp(t, 2)) {
	  if (0.90 <= gp(t, 2)) {
	    output(t) = "1/1";	  	  
	  } else {
	    output(t) = "./.";
	  }
	} else {
	  if (0.90 <= gp(t, 0)) {
	    output(t) = "0/0";
	  } else {
	    output(t) = "./.";
	  }
	}
      }
      // now, add on gp
      sprintf(buffer, ":%.3f,%.3f,%.3f:%.3f", gp(t, 0), gp(t, 1), gp(t, 2), gp(t, 1) + 2 * gp(t, 2));
      output(t) += buffer;
      if (use_read_proportions == 1) {
	sprintf(buffer, ":%.3f,%.3f,%.3f,%.3f",
		read_proportions(t, 0),
		read_proportions(t, 1),
		read_proportions(t, 2),
		read_proportions(t, 3));
	output(t) += buffer;	
      }
    }
    return output;
}
