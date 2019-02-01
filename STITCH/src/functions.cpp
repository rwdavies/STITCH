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



// EXAMPLE of how to export if desired
////' Reformat the reads for a single individual
////' 
////' @param x Describe parameters here.
////' @export
//// [[Rcpp::export]]






//' @export
// [[Rcpp::export]]
Rcpp::List rcpp_make_sampleReads_from_hap(const Rcpp::IntegerVector non_NA_cols, const int reference_phred, const Rcpp::IntegerVector reference_hap) {
    Rcpp::List sampleReads(non_NA_cols.length());
    for(int i = 0; i < non_NA_cols.length(); i++) {
        sampleReads[i]=Rcpp::List::create(0, non_NA_cols[i] - 1, reference_phred * (2 * reference_hap[i] - 1), non_NA_cols[i] - 1);
    }
    return sampleReads;
}




//' @export
// [[Rcpp::export]]
Rcpp::NumericVector increment2N(int yT, int xT, Rcpp::NumericVector y, Rcpp::NumericVector z) {
  Rcpp::NumericVector x(xT+1);
  int t;
  for(t=0; t<=yT-1; t++)
    x[z[t]]=x[z[t]]+y[t];
  return(x);
}



// why does this still exist?


//' @export
// [[Rcpp::export]]
Rcpp::List ram_test(
    const arma::mat& mat1,
    const Rcpp::NumericMatrix& mat2,
    arma::mat mat3,
    Rcpp::NumericMatrix mat4
) {
    double d1 = arma::accu(mat1);
    double d2 = Rcpp::sum(mat2);
    mat3(0, 0) = 3;
    mat4(0, 0) = 2;
    double d3 = arma::accu(mat3);
    double d4 = Rcpp::sum(mat4);
    arma::cube jimmy = arma::zeros(100, 1000, 2);
    return List::create(
                        Rcpp::Named("d1") = d1,
                        Rcpp::Named("d2") = d2,
                        Rcpp::Named("d3") = d3,
                        Rcpp::Named("d4") = d4
                        );
}





//' @export
// [[Rcpp::export]]
Rcpp::NumericVector get_random_values(int N) {
    Rcpp::NumericVector out = Rcpp::NumericVector(N);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    int n;
    double rand_uniform;
    for(n=0; n<N; n++){
        rand_uniform = dis(gen);        
        out(n)=rand_uniform;
    }
    return(out);
}




double print_times(
    double prev,
    int suppressOutput,
    std::string past_text,
    std::string next_text
) {
    if( suppressOutput == 0 ) {    
        double cur=clock();
        std::cout << std::setw (40) << past_text;
        printf ("- %.6f cpu sec -", ((double)cur - (double)prev)* 1.0e-6);
        std::cout << next_text << std::endl;
        prev=cur;
    }
    return prev;
}










//' @export
// [[Rcpp::export]]
Rcpp::List rcpp_make_fb_snp_offsets(
    const arma::mat& alphaHat_t,
    const arma::mat& betaHat_t,
    const arma::mat& blocks_for_output
) {
    int s, e;
    arma::mat alphaHatBlocks_t = arma::zeros(alphaHat_t.n_rows, blocks_for_output.n_rows);
    arma::mat betaHatBlocks_t = arma::zeros(betaHat_t.n_rows, blocks_for_output.n_rows);    
    for(int i_output=0; i_output < blocks_for_output.n_rows; i_output++) {
        s = blocks_for_output(i_output, 2); // these are 0-based. these are the grid entries
        e = blocks_for_output(i_output, 3);
        alphaHatBlocks_t.col(i_output) = alphaHat_t.col(s);
        betaHatBlocks_t.col(i_output) = betaHat_t.col(e);
    }
    return(wrap(Rcpp::List::create(
                                   Rcpp::Named("alphaHatBlocks_t") = alphaHatBlocks_t,
                                   Rcpp::Named("betaHatBlocks_t") = betaHatBlocks_t
                                   )));
}




//' @export
// [[Rcpp::export]]
arma::mat rcpp_calculate_fbd_dosage(
    const arma::mat& eHapsCurrent_t,
    const arma::mat& gamma_t,
    const Rcpp::IntegerVector& grid,
    const int snp_start_1_based,
    const int snp_end_1_based,
    const int grid_offset = 0
) {
    // basically, copy and paste, either using grid or not
    const int nSNPs = snp_end_1_based - snp_start_1_based + 1;
    const int K = eHapsCurrent_t.n_rows;
    arma::mat genProbs_t = arma::zeros(3, nSNPs);
    //arma::vec dosage = arma::zeros(nSNPs);
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
        eHapsCurrent_t_col = eHapsCurrent_t.col(t);
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
        genProbs_t(0, iSNP) = g0;
        genProbs_t(1, iSNP) = g1;
        genProbs_t(2, iSNP) = g2;        
        //for(j=1;j<=2;j++) {
        //    dosage(i_t) = dosage(i_t) + j * genProbs_t(j, i_t);
        //}
    }
    return(genProbs_t);
}











//' @export
// [[Rcpp::export]]
List cpp_read_reassign(
    arma::ivec ord,
    arma::ivec qnameInteger_ord,
    Rcpp::List sampleReadsRaw,
    int verbose,
    arma::ivec readStart_ord,
    arma::ivec readEnd_ord,
    int iSizeUpperLimit
) {
    // ord is 0-based original ordering
    // qnameInteger_ord is (ordered) integer representing reads
  
    // initialize
    int curRead = qnameInteger_ord[0];
    int iReadStart = 0;
    int nRawReads = sampleReadsRaw.size();
    int maxnSNPInRead = 1000;
    std::vector<int> base_bq(maxnSNPInRead);
    std::vector<int> base_pos(maxnSNPInRead); // there shouldnt be this many SNPs
    Rcpp::IntegerVector save_read(nRawReads); // over-sized
    Rcpp::LogicalVector save_this_read_check(nRawReads); // over-sized
    save_this_read_check.fill(false);
    int count = 0;
    int nReadsToSave = 0;
    int nSNPsInRead = -1;
    bool save_this_read;
    int iRead = 0;
    int j, r;
    arma::ivec bqL;
    arma::ivec posL;
    int to_add, k;

    // first loop, define which reads to save
    for (iRead = 0; iRead < nRawReads; iRead++ ) {
        if (qnameInteger_ord[iRead + 1] != curRead) {
            nSNPsInRead = -1;
            save_this_read = true;      
            for(j = iReadStart; j <= iRead; j++) {
                if (j < iRead) {
                    if ((readEnd_ord[j + 1] - readStart_ord[j]) > iSizeUpperLimit) {
                        save_this_read = false;
                    }
                }
            }
            if (save_this_read) {
                save_read(count) = curRead;
                save_this_read_check(iRead) = true;
                count++;
                nReadsToSave++;
            }
            iReadStart = iRead + 1;
            curRead = qnameInteger_ord[iRead + 1]; // + 1
        }
    }

    // second loop, use pre-defined sampleReads
    Rcpp::List sampleReads(nReadsToSave);
    
    curRead = qnameInteger_ord[0];    
    count = 0;
    iReadStart = 0;
    for (iRead = 0; iRead < nRawReads; iRead++ ) {
        if (qnameInteger_ord[iRead + 1] != curRead) {
            nSNPsInRead = -1;
            for(j = iReadStart; j <= iRead; j++) {
                r = ord[j];
                // so say first read is 0-based 0:2
                // want to take reads from ord[0:2]
                Rcpp::List readData = as<Rcpp::List>(sampleReadsRaw[r]);
                arma::ivec bqU = as<arma::ivec>(readData[2]);
                arma::ivec pRU = as<arma::ivec>(readData[3]);
                to_add = int(bqU.size());
                while ((nSNPsInRead + to_add) > maxnSNPInRead) {
                    base_bq.resize(maxnSNPInRead * 2);
                    base_pos.resize(maxnSNPInRead * 2);
                    maxnSNPInRead *= 2;
                }
                for(k = 0; k < to_add; k++) {
                    nSNPsInRead++;
                    base_bq[nSNPsInRead] = bqU[k];
                    base_pos[nSNPsInRead] = pRU[k];	  
                }
            }
            //
            if (save_this_read_check(iRead)) {
                arma::ivec bqL(nSNPsInRead + 1);
                arma::ivec posL(nSNPsInRead + 1);
                for(k = 0; k <= nSNPsInRead; k++) {
                    bqL[k] = base_bq[k];
                    posL[k] = base_pos[k];                    
                }
                //bqL = base_bq.subvec(0, nSNPsInRead);
                //posL = base_pos.subvec(0, nSNPsInRead);
                sampleReads[count] = Rcpp::List::create(nSNPsInRead, 0, bqL, posL);
                count++;
            }
            iReadStart = iRead + 1;
            curRead = qnameInteger_ord[iRead + 1]; // + 1
        }
    }
    
    Rcpp::List to_return = Rcpp::List::create(
        Rcpp::Named("sampleReads") = sampleReads,
        Rcpp::Named("save_read") = save_read
    );
    return to_return;
}





//' @export
// [[Rcpp::export]]
arma::mat make_gammaEK_t_from_gammaK_t(
    const arma::mat& gammaK_t,
    const int K,
    const Rcpp::IntegerVector& grid,
    const int snp_start_1_based,
    const int snp_end_1_based,
    const int grid_offset = 0
) {
    const int nSNPs = snp_end_1_based - snp_start_1_based + 1;
    arma::mat gammaEK_t = arma::zeros(K, nSNPs);
    int t;
    int prev_grid = -1;
    int iSNP;
    int cur_grid;
    arma::colvec gamma_t_col;
    for(iSNP = 0; iSNP < nSNPs; iSNP++) {
        t = iSNP + snp_start_1_based - 1;
        cur_grid = grid(t) - grid_offset;
        if (cur_grid > prev_grid) {
            gamma_t_col = gammaK_t.col(cur_grid);
            prev_grid = cur_grid;
        }
        gammaEK_t.col(iSNP) = gamma_t_col;
    }
    return(gammaEK_t);
}

