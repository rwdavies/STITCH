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


//' @export
// [[Rcpp::export]]
arma::rowvec test_eHaps_options(
    const arma::cube& cube_eHaps_t,
    const Rcpp::List& list_of_eHaps_t,
    const arma::mat& gamma_t,
    const arma::mat& eHaps_input,
    const std::string option,
    const int nSNPs,
    const int K,
    const int S
) {
    // true means cube
    arma::rowvec dosage = arma::zeros(1, nSNPs);
    int iSNP, k, s;
    arma::mat eHaps_t;
    if (option == "direct") {
      // simplest - work directly
        for(s = 0; s < S; s++) {
            for(iSNP = 0; iSNP < nSNPs; iSNP++) {
                for(k=0; k < K; k++) {
                    dosage(iSNP) += gamma_t(k, iSNP) * cube_eHaps_t(k, iSNP, s);
                }
            }
        }
    } else if (option == "direct_slice") {
      // simplest - work directly
        for(s = 0; s < S; s++) {
            for(iSNP = 0; iSNP < nSNPs; iSNP++) {
	      dosage(iSNP) += arma::sum(gamma_t.col(iSNP) % cube_eHaps_t.slice(s).col(iSNP));
            }
        }
	
    } else if (option == "list") {
        for(s = 0; s < S; s++) {
            eHaps_t = as<arma::mat>(list_of_eHaps_t[s]);
            for(iSNP = 0; iSNP < nSNPs; iSNP++) {
                for(k=0; k < K; k++) {
                    dosage(iSNP) += gamma_t(k, iSNP) * eHaps_t(k, iSNP);
                }
            }
        }
    } else if (option == "pre_get_slice") {
        for(s = 0; s < S; s++) {
            eHaps_t = cube_eHaps_t.slice(s);
            for(iSNP = 0; iSNP < nSNPs; iSNP++) {
                for(k=0; k < K; k++) {
                    dosage(iSNP) += gamma_t(k, iSNP) * eHaps_t(k, iSNP);
                }
            }
        }
    } else if (option == "best_possible") {
        for(s = 0; s < S; s++) {
            for(iSNP = 0; iSNP < nSNPs; iSNP++) {
	      dosage(iSNP) += arma::sum(gamma_t.col(iSNP) % eHaps_input.col(iSNP));
            }
        }
    }
    return(dosage);
}







