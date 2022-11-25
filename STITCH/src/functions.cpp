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
void get_min_from_position(
    const Rcpp::NumericVector & match_vec,
    const Rcpp::NumericVector & readStart_all,
    const Rcpp::NumericVector & readEnd_all,
    Rcpp::NumericVector & readMin,
    Rcpp::NumericVector & readMax
) {
    int w;
    for(int i=0; i<match_vec.length(); i++) {
        if (0 <= match_vec(i)) {
            w = match_vec(i);
            // if empty, take first value
            if (readMin(w) == -1) {
              readMin(w) = readStart_all(i);
            } else if (readStart_all(i) < readMin(w)) {
                readMin(w) = readStart_all(i);
            }
	    //                                    
            if (readMax(w) == -1) {
              readMax(w) = readEnd_all(i);
            } else if (readMax(w) < readEnd_all(i)) {
                readMax(w) = readEnd_all(i);
            }
        }
    }
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
Rcpp::LogicalVector rcpp_evaluate_bxtag(
    const arma::ivec& qnameInteger_ord,
    Rcpp::CharacterVector& bxtag_ord
) {
    const int n = bxtag_ord.length();
    Rcpp::LogicalVector bxtag_bad_ord = Rcpp::LogicalVector(n);
    bxtag_bad_ord.fill(false);
    bool perform_check;
    int iStart = 0;
    std::string doub_zero = "00";
    bool blank_local_bxtag;
    std::string first_bxtag, cur_bxtag;
    int n_non_doub_zero;
    int cur_qnameInteger = qnameInteger_ord(0);
    for(int iCur = 0; iCur < n; iCur++) {
        //std::cout << "iCur = " << iCur << " / " << n << std::endl;
        perform_check = false;
        if (iCur == (n - 1)) {
            perform_check = true;
        } else {
            if (qnameInteger_ord(iCur + 1) != cur_qnameInteger) {
                perform_check = true;
            }
        }
        //std::cout << "perform_check = " << perform_check << std::endl;        
        if (perform_check) {
            //std::cout << "------------------" << std::endl;
            first_bxtag = "";
            blank_local_bxtag = false;
            n_non_doub_zero = 0;
            for(int j = iStart; j <= iCur; j++) {
                cur_bxtag = bxtag_ord(j);
                //std::cout << "j = " << j << std::endl;
                //std::cout << "first_bxtag = " << first_bxtag << std::endl;                
                //std::cout << "cur_bxtag = " << cur_bxtag << std::endl;
                if (cur_bxtag.find(doub_zero) == std::string::npos) {
                    // not found i.e. it has no double 00 so can use
                    n_non_doub_zero++;
                    if (first_bxtag == "") {
                        first_bxtag = cur_bxtag;
                    } else {
                        if (cur_bxtag != first_bxtag) {
                            blank_local_bxtag = true;
                        }
                    }
                }
                //std::cout << "blank_local_bxtag = " << blank_local_bxtag << std::endl;        
            }
            // now afterwards, either blank out, or make the same
            if (blank_local_bxtag) {
                for(int j = iStart; j <= iCur; j++) {
                    bxtag_ord(j) = "00";
                    bxtag_bad_ord(j) = true;
                }
            }
            // also, if only one is zero
            if ((n_non_doub_zero > 0) & ((iCur - iStart + 1) != n_non_doub_zero)) {
                //std::cout << "apply rescue!" << std::endl;
                for(int j = iStart; j <= iCur; j++) {
                    bxtag_ord(j) = first_bxtag;
                }
            }
            // if all are zero
            if (n_non_doub_zero == 0) {
                for(int j = iStart; j <= iCur; j++) {
                    bxtag_bad_ord(j) = true;
                }
            }
            // reset
            if (iCur < (n - 1)) {
                cur_qnameInteger = qnameInteger_ord(iCur + 1);
                iStart = iCur + 1;
            }
        }
    }
    return(bxtag_bad_ord);
}




//' @export
// [[Rcpp::export]]
List cpp_read_reassign(
    const arma::ivec& ord,
    const arma::ivec& qnameInteger_ord,
    const arma::ivec& bxtagInteger_ord,
    const Rcpp::LogicalVector& bxtag_bad_ord,
    const Rcpp::CharacterVector& qname,
    const Rcpp::CharacterVector& bxtag,
    const Rcpp::CharacterVector& strand,
    const Rcpp::List& sampleReadsRaw,
    const arma::ivec& readStart_ord,
    const arma::ivec& readEnd_ord,
    const arma::ivec& readStart,
    const arma::ivec& readEnd,
    int iSizeUpperLimit,
    int bxTagUpperLimit,
    bool use_bx_tag,
    bool save_sampleReadsInfo = false,
    int maxnSNPInRead = 1000
) {
    // ord is 0-based original ordering
    // qnameInteger_ord is (ordered) integer representing reads
  
    // initialize
    int curRead = qnameInteger_ord[0];
    int curbxtag = bxtagInteger_ord[0];
    int iReadStart = 0;
    int nRawReads = sampleReadsRaw.size();
    std::vector<int> base_bq(maxnSNPInRead);
    std::vector<int> base_pos(maxnSNPInRead); // there shouldnt be this many SNPs
    //Rcpp::IntegerVector save_read(nRawReads); // over-sized
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
    long int minReadStart, maxReadEnd;
    std::string strand_to_out, qname_store;

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
                // record for each "read" with same qname
                for(j = iReadStart; j <= iRead; j++) {                
                    save_this_read_check(j) = true;
                }
                count++;
                nReadsToSave++;
            }
            iReadStart = iRead + 1;
            curRead = qnameInteger_ord[iRead + 1]; // + 1
        }
    }

    // second loop, use pre-defined sampleReads
    Rcpp::List sampleReads(nReadsToSave);
    Rcpp::CharacterVector sri_qname(nReadsToSave);
    Rcpp::CharacterVector sri_bxtag(nReadsToSave);
    Rcpp::CharacterVector sri_strand(nReadsToSave);    
    Rcpp::IntegerVector sri_minReadStart(nReadsToSave);
    Rcpp::IntegerVector sri_maxReadEnd(nReadsToSave);
    bool change_qname, save_condition_met, bx_tag_not_ok_to_continue, bx_tag_distance_exceeded, change_bx_tag;

    curRead = qnameInteger_ord[0];
    curbxtag = bxtagInteger_ord[0];
    count = 0;
    iReadStart = 0;
    for (iRead = 0; iRead < nRawReads; iRead++ ) {
        // if ((206 <= iRead) & (iRead <= 212)) {
        // std::cout << "------------------------" << std::endl;
        // std::cout << "iRead = " << iRead << std::endl;
        // std::cout << "iReadStart = " << iReadStart << std::endl;        
        // std::cout << "curRead = " << curRead << std::endl;                        
        // std::cout << "qnameInteger_ord[iRead + 1] = " << qnameInteger_ord[iRead + 1] << std::endl;
        // }
        //
        change_qname = qnameInteger_ord[iRead + 1] != curRead;
        // if ((206 <= iRead) & (iRead <= 212)) {        
        //     std::cout << "change_qname = " << change_qname << std::endl;
        // }
        if (!use_bx_tag) {
            save_condition_met = change_qname;
        } else {
            // if using bx tag, STOP (i.e. save) if
            // changing qname AND
            //   moving to a different bx tag OR
            //   bx tag not acceptable OR
            //   bx tag distance unacceptable (if qnames change AND distance exceeded)
            bx_tag_distance_exceeded = (change_qname) &                     \
                ((readStart_ord[iRead + 1] - readEnd_ord[iRead]) > bxTagUpperLimit);
            change_bx_tag = bxtagInteger_ord[iRead + 1] != curbxtag;
            bx_tag_not_ok_to_continue = (change_bx_tag) | \
                (bxtag_bad_ord[iRead]) | (bx_tag_distance_exceeded);
            save_condition_met = change_qname & bx_tag_not_ok_to_continue;
            //std::cout << "                      iRead = " << iRead << std::endl;
            //std::cout << "bxtagInteger_ord[iRead + 1] = " << bxtagInteger_ord[iRead + 1] << std::endl;            
            //std::cout << "                   curbxtag = " << curbxtag << std::endl;                        
            // std::cout << "change_bx_tag = " << change_bx_tag << std::endl;
            // std::cout << "bxtag_bad_ord[iRead] = " << bxtag_bad_ord[iRead] << std::endl;
            // std::cout << "bx_tag_distance_exceeded = " << bx_tag_distance_exceeded << std::endl;
            //std::cout << "  bx_tag_not_ok_to_continue = " << bx_tag_not_ok_to_continue << std::endl;            
            // std::cout << "change_qname = " << change_qname << std::endl;
            //std::cout << "         save_condition_met = " << save_condition_met << std::endl;
        }
        if (change_qname) {
            curRead = qnameInteger_ord[iRead + 1]; // + 1
        }
        if (save_condition_met) {
            //std::cout << "SAVING WITH iReadStart = " << iReadStart << ", iRead = " << iRead << std::endl;
            nSNPsInRead = -1; // this is 0-based (why did I do this)
            minReadStart = 2147483647;
            maxReadEnd = -1;
            strand_to_out = "";
            if (use_bx_tag) {
                qname_store="";
            }
            for(j = iReadStart; j <= iRead; j++) {
                // check if we should be saving this bit
                //std::cout << "j = " << j << ", save_this_read_check(j) = " << save_this_read_check(j) << std::endl;
                if (save_this_read_check(j)) {
                    r = ord[j];
                    if (use_bx_tag) {
                        if (qname_store == "") {
                            qname_store = qname[r];
                        } else {
                            qname_store = qname_store + "," + qname[r];
                        }
                    }
                    // so say first read is 0-based 0:2
                    // want to take reads from ord[0:2]
                    Rcpp::List readData = as<Rcpp::List>(sampleReadsRaw[r]);
                    arma::ivec bqU = as<arma::ivec>(readData[2]);
                    arma::ivec pRU = as<arma::ivec>(readData[3]);
                    to_add = int(bqU.size());
                    while ((nSNPsInRead + 1 + to_add) > maxnSNPInRead) {
                        base_bq.resize(maxnSNPInRead * 2);
                        base_pos.resize(maxnSNPInRead * 2);
                        maxnSNPInRead *= 2;
                    }
                    for(k = 0; k < to_add; k++) {
                        nSNPsInRead++;
                        base_bq[nSNPsInRead] = bqU[k];
                        base_pos[nSNPsInRead] = pRU[k];	  
                    }
                    // can also take from r otherwise                
                    if (save_sampleReadsInfo) {
                        if (readStart[r] < minReadStart) {
                            minReadStart = readStart[r];
                        }
                        if (maxReadEnd < readEnd[r]) {
                            maxReadEnd = readEnd[r];
                        }
                        strand_to_out.append(strand[r]);
                    }
                }
            }
            //
            if (nSNPsInRead >= 0) {
                arma::ivec bqL(nSNPsInRead + 1);
                arma::ivec posL(nSNPsInRead + 1);
                for(k = 0; k <= nSNPsInRead; k++) {
                    bqL[k] = base_bq[k];
                    posL[k] = base_pos[k];
                }
                //bqL = base_bq.subvec(0, nSNPsInRead);
                //posL = base_pos.subvec(0, nSNPsInRead);
                sampleReads[count] = Rcpp::List::create(nSNPsInRead, 0, bqL, posL);
                // can draw from ord(iRead);
                if (save_sampleReadsInfo) {
                    if (!use_bx_tag) {
                        sri_qname[count] = qname[r];
                    } else {
                        sri_qname[count] = qname_store;
                        sri_bxtag[count] = bxtag[r];
                    }
                    sri_minReadStart[count] = minReadStart;
                    sri_maxReadEnd[count] = maxReadEnd;
                    sri_strand[count] = strand_to_out;
                }
                // if ((206 <= iRead) & (iRead <= 212)) {
                //     std::cout << " saving bit" << std::endl;
                //     std::cout << "count = " << count << std::endl;                    
                // }
                count++;
            }
            iReadStart = iRead + 1;
            curbxtag = bxtagInteger_ord[iRead + 1]; // + 1
        }
    }

    Rcpp::DataFrame sampleReadsInfo;
    if (save_sampleReadsInfo) {
        if (!use_bx_tag) {
            sampleReadsInfo = Rcpp::DataFrame::create(
                Rcpp::Named("qname") = sri_qname,
                Rcpp::Named("strand") = sri_strand,
                Rcpp::Named("minReadStart") = sri_minReadStart,
                Rcpp::Named("maxReadEnd") = sri_maxReadEnd
            );
        } else {
            sampleReadsInfo = Rcpp::DataFrame::create(
                Rcpp::Named("qname") = sri_qname,
                Rcpp::Named("bxtag") = sri_bxtag,                
                Rcpp::Named("strand") = sri_strand,
                Rcpp::Named("minReadStart") = sri_minReadStart,
                Rcpp::Named("maxReadEnd") = sri_maxReadEnd
            );
        }
    }
    Rcpp::List to_return = Rcpp::List::create(
        Rcpp::Named("sampleReads") = sampleReads,
        Rcpp::Named("sampleReadsInfo") = sampleReadsInfo
    );
    //         Rcpp::Named("save_read") = save_read,
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
    double& prev,
    const int suppressOutput,
    std::string& prev_section,
    std::string& next_section,
    const int grid_offset = 0
) {
    next_section="make gammaEK_t from gammaK_t";
    prev=print_times(prev, suppressOutput, prev_section, next_section);
    prev_section=next_section;
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

