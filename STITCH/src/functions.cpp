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




////Rcpp::List pseudoHaploid_update_model_9(const arma::vec& pRgivenH1, const arma::vec& pRgivenH2, const arma::mat& eMatHap_t1, const arma::mat& eMatHap_t2, const arma::mat& gamma_t1, const arma::mat& gamma_t2, int K, const arma::ivec& srp) {

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
Rcpp::List reformatReads(arma::ivec mapq, arma::ivec readStart, arma::ivec readEnd, const int numberOfReads, const int T, arma::ivec L, Rcpp::CharacterVector seqRead, Rcpp::CharacterVector qualRead,  Rcpp::CharacterVector ref,  Rcpp::CharacterVector alt, Rcpp::List splitCigarRead, arma::ivec lengthOfSplitCigarRead, const int bqFilter, const int verbose) {
//
// new variables
//
  int x1, x2, y, iM, iRead, t;
  int nSNPInRead = -1;
  int nReadSpanningSNPs = 0;
  char s;
  Rcpp::List sampleReads;
  arma::ivec seqLocal(1000);
  arma::ivec qualLocal(1000);
  arma::ivec posLocal(1000); // there shouldnt be this many SNPs
  int refPosition, refOffset, strandOffset;
  int iNumOfMs, cigarLength;
  Rcpp::CharacterVector cigarType(1);
  Rcpp::CharacterVector cigar2(1);
  int whileVar, localbq;
  //
  // begin
  //
  // want - tMin - first SNP after last SNP before read
  // also - tMax - first SNP before last SNP after read
  // if the distance is >=0 between them, it is in that read
  //
  int tMin=0; 
  int tMax=-1; // left and right boundaries of what SNPs to look at
  for(iRead=0; iRead<=numberOfReads-1; iRead++)
  {
    // get current read position
    whileVar=0;
    // for t-Min
    while(whileVar==0)
    {
      // dont continue if too far
      if(tMin<(T-1))
      {
        // continue 
        if(L(tMin)<readStart[iRead])
        {
          tMin++;
        } else {
          whileVar=1; // break loop - done!
        }
      } else {
        whileVar=1; // break loop
      }
    }
    // now - go right
    whileVar=0;
    while(whileVar==0)
    {
      // dont continue if too far
      if(tMax < (T - 1))
      {
	if(L(tMax + 1) <= readEnd[iRead])
        {
	  // number is only to deal with weird cigars	  
          tMax++;
        } else {
          whileVar=1; // break loop - done!
        }
      } else {
	tMax = T - 1; // done, cannot be larger
        whileVar=1; // break loop
      }
    }
    //
    //
    // now, for this read, calculate whether there are SNPs - only bother if tMin>=tMax
    //
    //
    if(verbose == 1) {
            std::cout << "iRead=" << iRead;
            std::cout << ",tMin=" << tMin;
            std::cout << ",tMax=" << tMax;
            std::cout << "\\n";
    }
    nSNPInRead=-1;
    if(tMin<=tMax)
    {
      refPosition=readStart(iRead);
      Rcpp::List cigarReadInfo = as<Rcpp::List>(splitCigarRead(iRead));
      iNumOfMs = lengthOfSplitCigarRead(iRead);
      // set some things
      refOffset=0; // offset against the reference sequence
      strandOffset=0; // offset in the strand
      // get cigar info from the read
      arma::ivec cigarLengthVec = as<arma::ivec>(cigarReadInfo(0));
      Rcpp::CharacterVector cigarTypeVec = as<Rcpp::CharacterVector>(cigarReadInfo(1));
      //
      // now, loop over each part of the read (M, D=del, I=ins)
      //
      for(iM=0;iM<=iNumOfMs;iM++)
      {
        cigarLength=cigarLengthVec(iM);
        cigarType(0)=cigarTypeVec(iM);
        // if its an M - scan
        if(cigarType(0)=="M")
        {
          x1 = refPosition + refOffset; // left part of M
          x2 = refPosition + refOffset + cigarLength-1; // right part of M
          for(t=tMin; t<=tMax; t++) // determine whether that snps is spanned by the read
          {
            y = L[t];
            if(x1 <= y && y <= x2) // if this is true - have a SNP!
            {
              s = seqRead[iRead][y-refPosition-refOffset+strandOffset];
              // check if ref or ALT - only keep if true
              // also only use if BQ at least bqFilter (17) (as in 17 or greater)
              localbq=int(qualRead[iRead][y-refPosition-refOffset+strandOffset])-33;
              // also bound BQ above by MQ
              if(localbq>mapq(iRead)) // if greater, than reduce
                localbq=mapq(iRead);
//    if(verbose == 1) {
//            std::cout << "checking - iRead=" << iRead;
//            std::cout << ",s=" << s;
//            std::cout << ",ref[t][0]=" << ref[t][0];
//            std::cout << ",alt[t][0]=" << alt[t][0];
//            std::cout << ",t=" << t;
//           std::cout << ",y=" << y;
//            std::cout << ",x1=" << x1;
//            std::cout << ",x2=" << x2;
//            std::cout << "\\n";
//    }
              if((s==ref[t][0] || s==alt[t][0]) && (localbq>=bqFilter))
              {
                // is this the reference or alternate?
                nSNPInRead = nSNPInRead+1;
                if(s==ref[t][0])
                {
                  seqLocal[nSNPInRead] = 0;
                  qualLocal[nSNPInRead] = - localbq;
                }
                if(s==alt[t][0])
                {
                  seqLocal[nSNPInRead] = 1;
                  qualLocal[nSNPInRead] = localbq;
                }
                posLocal[nSNPInRead] = t;
              } // end of check if ref or alt
            } // end of whether this SNP intersects read
          } // end of loop on SNP
          // now, bump up ref and pos offset by x1
          refOffset=refOffset + cigarLength;
          strandOffset=strandOffset + cigarLength;
        } // end of if statement on whether cigar type is M
        // if it is an insertion - bump the strand offset
        if(cigarType(0)=="I")
          strandOffset=strandOffset+cigarLength;
        // if it is a deletion - bump the reference position
        if(cigarType(0)=="D")
          refOffset=refOffset+cigarLength;
      } // close loop on M
    } // end of check on whether there can be results to run
    //
    // save result if it is worth saving!
    //
    if(nSNPInRead > -1) // save result!
    {
    //if(verbose == 1) {
    //        std::cout << "saving - iRead=" << iRead;
     //       std::cout << ",nSNPInRead=" << nSNPInRead;
     //       std::cout << ",tMin=" << tMin;
     //       std::cout << ",tMax=" << tMax;
     //       std::cout << "\\n";
   // }
      arma::ivec pR = posLocal.subvec(0,nSNPInRead); // position (R means Read)
      arma::ivec sR = seqLocal.subvec(0,nSNPInRead); // sequence
      arma::ivec qR = qualLocal.subvec(0,nSNPInRead); // quality
      //
      // save results but dont label list elements to save space
      //
      // save a smaller version unless need to debug
      sampleReads.push_back(Rcpp::List::create(pR.size()-1,0,qR,pR,iRead));
      nReadSpanningSNPs = nReadSpanningSNPs +1;
    }  // end of save result
  } // end of loop on read
  //
  // done
  //
  return(wrap(sampleReads));
}












double print_times (double prev, int suppressOutput, std::string past_text, std::string next_text) {
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
Rcpp::List forwardBackwardDiploid(const Rcpp::List& sampleReads,const int nReads, const arma::vec& pi, const arma::mat& transMatRate_t, const arma::mat& alphaMat_t, const arma::mat& eHaps_t, const double maxDifferenceBetweenReads, const int whatToReturn, const int Jmax, const int suppressOutput) {
  double prev=clock();
  std::string prev_section="Null";
  std::string next_section="Initialize variables";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  //
  // constants
  //
  const int T = eHaps_t.n_cols;
  const int K = eHaps_t.n_rows; // traditional K for haplotypes
  const int KK = K*K; // KK is number of states / traditional K for HMMs
  //
  // new variables
  //
  // variables working on nReads
  arma::mat alphaHat_t = arma::zeros(KK,T);  
  arma::mat betaHat_t = arma::zeros(KK,T);
  arma::rowvec c = arma::zeros(1,T);
  arma::mat gamma_t = arma::zeros(KK,T);  
  arma::mat eMat_t = arma::ones(KK,T);
  // define eMatHap to work on the number of reads
  arma::mat eMatHap_t = arma::ones(K,nReads);  
  // variables for faster forward backward calculation
  double alphaConst, betaConst;
  arma::vec alphaTemp1 = arma::zeros(K);
  arma::vec alphaTemp2 = arma::zeros(K);
  arma::vec betaTemp1 = arma::zeros(K);
  arma::vec betaTemp2 = arma::zeros(K);
  // variables working on full space
  arma::mat jUpdate_t = arma::zeros(K,T-1);  
  arma::cube gammaUpdate_t = arma::zeros(K,T,2);
  arma::mat hapProbs = arma::zeros(T,K);
  arma::mat genProbs = arma::zeros(T,3);
  arma::vec dosage = arma::zeros(T);
  int jj, j, k, kk, t, k1, k2, iRead;
  int J, readSNP;
  double kl1, kl2;        
  double x, a, b, d, e, d1, d2, d3, eps;
  double pR = 0;
  double pA = 0;
  //
  //
  // eMat - make complete
  //
  //
  next_section="Initialize eMatHap";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  //
  for(iRead=0; iRead<=nReads-1; iRead++)
  {
    Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
    // recal that below is what is used to set each element of sampleRead
    // note - this is no longer quite accurate
    // sampleReads.push_back(Rcpp::List::create(nU,d,phiU,pRU));
    J = as<int>(readData[0]); // number of Unique SNPs on read
    readSNP = as<int>(readData[1]); // leading SNP from read
    arma::ivec bqU = as<arma::ivec>(readData[2]); // bq for SNPs
    arma::ivec pRU = as<arma::ivec>(readData[3]); // position of each SNP from 0 to T-1
    // once each SNP is done, have P(read | k), can multiply to get P(read|(k1,k2))
    if(J>=Jmax)
      J=Jmax;
    // for each SNP in the read
    for(j=0; j<=J; j++) 
    {
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
      jj=pRU(j);
      for(k=0; k<=K-1; k++)
          eMatHap_t(k,iRead) = eMatHap_t(k,iRead) * \
              ( eHaps_t(k,jj) * pA + (1 - eHaps_t(k,jj)) * pR);
    }
    //
    // cap P(read|k) to be within maxDifferenceBetweenReads orders of magnitude
    //
    x=0;
    for(k=0; k<=K-1; k++)
        if(eMatHap_t(k,iRead)>x)
            x=eMatHap_t(k,iRead);
    // x is the maximum now
    for(k=0; k<=K-1; k++)
        if(eMatHap_t(k,iRead)<x/maxDifferenceBetweenReads)
            eMatHap_t(k,iRead)=x/maxDifferenceBetweenReads;
  }
  //
  // once we have all the eMatHaps, ie probabilities from reads, make eMat from this
  //
  next_section="Initialize eMat";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  d=maxDifferenceBetweenReads*maxDifferenceBetweenReads;
  //
  for(iRead=0; iRead<=nReads-1; iRead++) {
      Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
      readSNP = as<int>(readData[1]); // leading SNP from read
      for(k1=0; k1<=K-1; k1++) {
          for(k2=0; k2<=K-1; k2++) {
              // multiply by former value if >=2 reads at a locus
              eMat_t(k1+K*k2,readSNP) = eMat_t(k1+K*k2,readSNP) * (0.5 * eMatHap_t(k1,iRead) + 0.5 * eMatHap_t(k2,iRead));
          }
      }  // end of SNP in read
  }
  //
  // cap eMat, ie P(reads | k1,k2) to be within maxDifferenceBetweenReads^2
  //
  for(iRead=0; iRead<=nReads-1; iRead++) {
    Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
    readSNP = as<int>(readData[1]); // leading SNP from read
    // first, get maximum
    x=0;
    for(k=0; k<=KK-1; k++)
        if(eMat_t(k,readSNP)>x)
            x=eMat_t(k,readSNP);
    // x is the maximum now
    for(k=0; k<=KK-1; k++)
        if(eMat_t(k,readSNP)<(x/d))
            eMat_t(k,readSNP)=x/d;
  } // end of loop on t
  //
  //
  // forward recursion
  //
  //
  next_section="Forward recursion";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  for(k1=0; k1<=K-1; k1++)
    for(k2=0; k2<=K-1; k2++)
        alphaHat_t(k1+K*k2,0) = pi(k1) * pi(k2) * eMat_t(k1+K*k2,0);
  c(0) = 1 / sum(alphaHat_t.col(0));
  alphaHat_t.col(0) = alphaHat_t.col(0) * c(0);
  //
  //
  // forward recursion
  //
  //
  for(t=1; t<=T-1; t++) {
    // calculate necessary things
    alphaTemp1.fill(0);
    alphaTemp2.fill(0);
    for(k1=0; k1<=K-1; k1++) {
      for(k2=0; k2<=K-1; k2++) {
          kk=k1+K*k2;
          alphaTemp1(k2) = alphaTemp1(k2) + alphaHat_t(kk,t-1);
          alphaTemp2(k1) = alphaTemp2(k1) + alphaHat_t(kk,t-1);
      }
    }
    // now make constant over whole thing
    alphaConst=0;
    for(kk=0; kk<=KK-1; kk++)
        alphaConst = alphaConst + alphaHat_t(kk,t-1);
    alphaConst = alphaConst * transMatRate_t(2,t-1);
    for(k=0; k<=K-1; k++) {
        alphaTemp1(k)=alphaTemp1(k) * transMatRate_t(1, t-1);
        alphaTemp2(k)=alphaTemp2(k) * transMatRate_t(1, t-1);
    }
    // 
    for(k1=0; k1<=K-1; k1++) {
        for(k2=0; k2<=K-1; k2++) {
            kk=k1+K*k2;
            alphaHat_t(kk,t) = eMat_t(kk,t) *    \
                (alphaHat_t(kk,t-1) * transMatRate_t(0,t-1) + \
                 (alphaMat_t(k1,t-1) * alphaTemp1(k2) +    \
                  alphaMat_t(k2,t-1) * alphaTemp2(k1)) + \
                 alphaMat_t(k1,t-1) * alphaMat_t(k2,t-1) * alphaConst);
        }
    }
    // do scaling now
    c(t) = 1 / sum(alphaHat_t.col(t));
    alphaHat_t.col(t) = alphaHat_t.col(t) * c(t);
  }
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
  betaHat_t.col(T-1).fill(c(T-1));
  for(t = T-2; t>=0; --t) {
    betaTemp1.fill(0);
    betaTemp2.fill(0);
    betaConst=0;
    for(k1=0; k1<=K-1; k1++) {
        for(k2=0; k2<=K-1; k2++) {
            kk=k1+K*k2;
            d=betaHat_t(kk,t+1) * eMat_t(kk,t+1);
            betaTemp1(k1) = betaTemp1(k1) + \
                d * alphaMat_t(k2,t);
            betaTemp2(k2) = betaTemp2(k2) + \
                d * alphaMat_t(k1,t);
            betaConst = betaConst + \
                d * alphaMat_t(k1,t) * alphaMat_t(k2,t);
        }
    }
    // add transMatRate to constants
    d=transMatRate_t(1,t);
    for(k=0; k<=K-1; k++) {
        betaTemp1(k) = betaTemp1(k) * d;
        betaTemp2(k) = betaTemp2(k) * d;
    }
    betaConst = betaConst * transMatRate_t(2,t);
    // final calculation
    for(k1=0; k1<=K-1; k1++) {
        for(k2=0; k2<=K-1; k2++) {
            kk=k1+K*k2;
            betaHat_t(kk,t) = eMat_t(kk,t+1) * betaHat_t(kk,t+1) *\
                transMatRate_t(0,t) +         \
                betaTemp1(k1) + betaTemp2(k2) + \
                betaConst;
        }
    }
    // apply scaling
    betaHat_t.col(t) = betaHat_t.col(t) * c(t);
  }
  //
  // DONE looping!
  //
  //
  // make gamma
  //
  next_section="Make gamma";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  //
  gamma_t = alphaHat_t % betaHat_t;
  //for(kk=0; kk <= KK - 1; kk++)
  //    for(t=0; t<= T-1; t++)
  //        gamma(t,kk) = alphaHat_t(kk,t) * betaHat_t(kk,t);
  //
  for(t=0; t<= T-1; t++)
      gamma_t.col(t) = gamma_t.col(t) / c(t);
  //
  //
  // get dosages
  //  - only needed if final run
  //
  //
  next_section="Get dosages";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  if((whatToReturn == 2) | (whatToReturn==0)) // final run or debugging only
  {
    for(t=0; t<=T-1; t++)
    {
      for(k1=0; k1<=K-1;k1++)
      {
        for(k2=0; k2<=K-1;k2++)
        {
          k=k1+K*k2;
          hapProbs(t,k1)= hapProbs(t,k1) + gamma_t(k,t)/2;
          hapProbs(t,k2)= hapProbs(t,k2) + gamma_t(k,t)/2;
          a=eHaps_t(k1, t);
          b=eHaps_t(k2, t);
          genProbs(t,0) = genProbs(t,0) + \
              gamma_t(k,t) * (1-a) * (1-b);
          genProbs(t,1) = genProbs(t,1) + \
              gamma_t(k,t) * (a * (1 - b) + (1 - a) * b);
          genProbs(t,2) = genProbs(t,2) + \
              gamma_t(k,t) * a * b;
        } // k2
      } //k1
      for(j=1;j<=2;j++)
        dosage(t)=dosage(t)+j*genProbs(t,j)/2;
    } // end of loop on t
  } // end of if statement
  //
  //
  // do gamma update here - save large matrix, just put in necessary values
  //
  //
  next_section="Gamma update";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  //
  // only, mathematically correct version
  //
  for(iRead=0; iRead<=nReads-1; iRead++)
  {
      // recal that below is what is used to set each element of sampleRead
      // sampleReads.push_back(Rcpp::List::create(nU,d,phiU,pRU));
      Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
      int J = as<int>(readData[0]); // number of SNPs on read
      int cr = as<int>(readData[1]); // central SNP
      arma::ivec bqU = as<arma::ivec>(readData[2]); // bq for each SNP
      arma::ivec pRU = as<arma::ivec>(readData[3]); // position of each SNP from 0 to T-1
      for(j=0; j<=J; j++) // for each SNP in the read
      {
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
        //
        // RECAL  eMatHap(iRead,k) = eMatHap(iRead,k) * ( eHaps(pRU(j),k) * pA + (1-eHaps(pRU(j),k)) * pR);
        // RECAL  eMat(readSNP,k1+K*k2) = eMat(readSNP,k1+K*k2) * (0.5 * eMatHap(iRead,k1) + 0.5 * eMatHap(iRead,k2));
        //
        //
        for(k=0; k<=K-1;k++) {
            d1 = pA * eHaps_t(k,t);
            d2 = pR * (1-eHaps_t(k,t));
            d3 = d1 / (d1 + d2); // this is all I need
            a = eMatHap_t(k,iRead);
            for(k1=0; k1<=K-1; k1++) {
                kl1 = k+K*k1;
                kl2 = k1+K*k;            
                b = eMatHap_t(k1, iRead);
                e = (gamma_t(kl1,cr) + gamma_t(kl2,cr)) * ( a / (a + b));
                gammaUpdate_t(k,t,0) = gammaUpdate_t(k,t,0) + e * d3;
                gammaUpdate_t(k,t,1) = gammaUpdate_t(k,t,1) + e;
            } // loop on other haplotype
        } // end of loop onk
      } // end of loop on SNP within read
  }
  //
  //
  // make jUpdate
  //
  //
  next_section="Make xi-like calculations";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  for(t=0; t<=T-2; t++) {
      alphaTemp1.fill(0);    
      for(k1=0; k1<=K-1; k1++) {
          for(k2=0; k2<=K-1; k2++) {
              alphaTemp1(k2) = alphaTemp1(k2) + alphaHat_t(k1+K*k2,t);
          }
      }
      //
      // now do proper calculation
      //
      for(k1=0; k1<=K-1; k1++) {
          for(k2=0; k2<=K-1; k2++) {
              kk=k1+K*k2;
              jUpdate_t(k1,t) = jUpdate_t(k1,t) +           \
                  (transMatRate_t(1,t) * alphaTemp1(k2) +  \
                   transMatRate_t(2,t) * alphaMat_t(k2,t)) * \
                  betaHat_t(kk,t+1) * eMat_t(kk,t+1);
          }
          jUpdate_t(k1,t) = jUpdate_t(k1,t) * 2 * alphaMat_t(k1,t);
      }   // end of loop on k
  } // end of loop on t
  //
  //
  // done - return necessary information
  //
  //
  next_section="Done";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  if(whatToReturn==0) // full, debugging return
    return(wrap(Rcpp::List::create(
      Rcpp::Named("hapProbs") = hapProbs,
      Rcpp::Named("dosage") = dosage,
      Rcpp::Named("genProbs") = genProbs,
      Rcpp::Named("gammaUpdate_t") = gammaUpdate_t,
      Rcpp::Named("gamma_t") = gamma_t,      
      Rcpp::Named("betaHat_t") = betaHat_t,
      Rcpp::Named("alphaHat_t") = alphaHat_t,      
      Rcpp::Named("jUpdate_t") = jUpdate_t,      
      Rcpp::Named("eMat_t") = eMat_t,
      Rcpp::Named("eMatHap_t") = eMatHap_t)));
  if(whatToReturn==1) // update
    return(wrap(Rcpp::List::create(
      Rcpp::Named("gammaUpdate_t") = gammaUpdate_t,
      Rcpp::Named("gamma_t") = gamma_t,
      Rcpp::Named("jUpdate_t") = jUpdate_t)));
  if(whatToReturn==2) // final run - no updating needed
    return(wrap(Rcpp::List::create(
      Rcpp::Named("hapProbs") = hapProbs,
      Rcpp::Named("gamma_t") = gamma_t,
      Rcpp::Named("dosage") = dosage,
      Rcpp::Named("genProbs") = genProbs,
      Rcpp::Named("jUpdate_t") = jUpdate_t,
      Rcpp::Named("gammaUpdate_t") = gammaUpdate_t)));
  // just put something to remove warning
  return(wrap(Rcpp::List::create(
      Rcpp::Named("c") = c)));
  
}















//' @export
// [[Rcpp::export]]
Rcpp::List forwardBackwardHaploid(const Rcpp::List& sampleReads, const int nReads, const arma::vec pi, const arma::mat& transMatRate_t, const arma::mat& alphaMat_t, const arma::mat& eHaps_t, const double maxDifferenceBetweenReads, const int whatToReturn, const int Jmax, const int suppressOutput, const int model, const arma::vec& pRgivenH1, const arma::vec& pRgivenH2, arma::mat pState) {
  double prev=clock();
  std::string prev_section="Null";
  std::string next_section="Initialize variables";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  //
  // constants
  //
  //
  const int T = eHaps_t.n_cols;
  const int K = eHaps_t.n_rows; // traditional K for haplotypes
  //
  // new variables
  //
  // variables working on nReads
  arma::mat alphaHat_t = arma::zeros(K,T);  
  arma::mat betaHat_t = arma::zeros(K,T);
  arma::rowvec c = arma::zeros(1,T);  
  arma::mat gamma_t = arma::zeros(K,T);
  // define eMatHap to work on the number of reads
  arma::mat eMatHapOri_t = arma::ones(K,nReads);
  arma::mat eMatHap_t = arma::zeros(K,nReads);
  // eMatHapSNP works on the SNPs themselves
  arma::mat eMatHapSNP_t = arma::ones(K,T);
  // variables working on full space
  arma::mat jUpdate_t = arma::zeros(K,T-1);
  arma::cube gammaUpdate_t = arma::zeros(K,T,2);
  // variables for transition matrix and initialization
  // int variables and such
  int j, jj, k, t, k1, iRead;
  int J, readSNP;
  double x, b, d, d1, d2, d3, a1, a2, y;
  double alphaConst;
  double eps;
  double pRead = 0;
  double pR = 0;
  double pA = 0;
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
  for(iRead=0; iRead<=nReads-1; iRead++)
  {
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
            eMatHapOri_t(k,iRead) = eMatHapOri_t(k,iRead) *   \
                ( eHaps_t(k, jj) * pA + (1-eHaps_t(k,jj)) * pR);
    } 
    // 
    //
    x=pRgivenH1(iRead) / (pRgivenH1(iRead) + pRgivenH2(iRead));
    //
    if((model==7) | (model==9))
        for(k=0; k<=K-1; k++)
            eMatHap_t(k,iRead) = x * eMatHapOri_t(k,iRead) + (1-x) * pRgivenH2(iRead);
    if((model==8) | (model==10)) {
        if(model==8) {
            pRead=0;
            for(k=0; k<=K-1; k++)
                pRead=pRead + eMatHapOri_t(k,iRead) * pState(readSNP,k);
        }
        for(k=0; k<=K-1; k++)
            eMatHap_t(k,iRead) = x * eMatHapOri_t(k,iRead) + (1-x) * pRead;
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
  //
  //
  // once we have eMatHap, ie probabilities from reads, make eMatHapSNPsfrom this
  //
  //
  next_section="Make eMat";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  //
  for(iRead=0; iRead<=nReads-1; iRead++) {
    Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
    readSNP = as<int>(readData[1]); // leading SNP from read
    for(k=0; k<=K-1; k++) {
        eMatHapSNP_t(k, readSNP) = eMatHapSNP_t(k, readSNP) *  \
            eMatHap_t(k,iRead);
    }
  }
  //
  //
  // afterwards - cap per-SNP by difference squared 
  //
  //
  d = maxDifferenceBetweenReads * maxDifferenceBetweenReads;
  // now - afterward - cap eMatHapSNP
  for(t=0; t<=T-1; t++) {
      // if eMatHapSNP(t, 0) != 0, proceed
      if (eMatHapSNP_t(0, t) > 0) {
          x=0;
          for(k=0; k<=K-1; k++)
              if(eMatHapSNP_t(k, t)>x)
                  x=eMatHapSNP_t(k, t);
          // x is the maximum now
          x=x/d;
          for(k=0; k<=K-1; k++)
              if(eMatHapSNP_t(k, t)<x)
                  eMatHapSNP_t(k, t)=x;
      }
  }
  //
  //
  // forward recursion
  //
  //
  next_section="Forward recursion";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
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
          alphaHat_t(k,t) =  eMatHapSNP_t(k,t) *            \
              ( transMatRate_t(0, t-1) * alphaHat_t(k,t-1) +   \
                alphaConst * alphaMat_t(k,t-1));
      //
      c(t) = 1 / sum(alphaHat_t.col(t));
      alphaHat_t.col(t) = alphaHat_t.col(t) * c(t);
  }
  //
  //
  // backward recursion
  //
  //
  next_section="Backward recursion";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  betaHat_t.col(T-1).fill(c(T-1));
  for(t = T-2; t>=0; --t) {
      x = 0;
      for(k=0; k<=K-1; k++)
          x = x + alphaMat_t(k,t) * eMatHapSNP_t(k,t+1) * betaHat_t(k,t+1);
      x = x * transMatRate_t(1, t);
      for(k=0; k<=K-1; k++)
          betaHat_t(k,t) = x + \
              transMatRate_t(0, t) * eMatHapSNP_t(k,t+1) * betaHat_t(k,t+1);
      // 
      betaHat_t.col(t) = betaHat_t.col(t) * c(t);
  }
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
  for(t=0; t<= T-1; t++)
      gamma_t.col(t) = gamma_t.col(t) / c(t);
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
  //
  // I should probably delete the below
  //
  if((model==7) | (model==8))
  {
    for(iRead=0; iRead<=nReads-1; iRead++)
    {
      // recal that below is what is used to set each element of sampleRead
      // sampleReads.push_back(Rcpp::List::create(nU,d,phiU,pRU));
      Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
      int J = as<int>(readData[0]); // number of SNPs on read
      arma::ivec bqU = as<arma::ivec>(readData[2]); // bq for each SNP
      arma::ivec pRU = as<arma::ivec>(readData[3]); // position of each SNP from 0 to T-1
      for(j=0; j<=J; j++) // for each SNP in the read
      {
        t=pRU(j); // position of this SNP in full T sized matrix
        //
        // first haplotype (ie (k,k1))
        //
        for(k=0; k<=K-1;k++)
        {
          // less than 0 - reference base more likely
          if(bqU(j)<0)
          {
            eps = pow(10,(double(bqU(j))/10));
            pR=1-eps;
            pA=eps/3;
          }
          if(bqU(j)>0)
          {
            eps = pow(10,(-double(bqU(j))/10));
            pR=eps/3;
            pA=1-eps;
          }
          //
          //
          //    
          x=pRgivenH1(iRead) / (pRgivenH1(iRead) + pRgivenH2(iRead));
          d1 = x * pA * eHaps_t(k, t);
          d2 = x * pR * (1-eHaps_t(k, t));
          d= gamma_t(k,t) / eMatHap_t(k,iRead) ;
          // d1 = x * phiU(j) * eHaps(t,k);
          // d2 = x * (1-phiU(j)) * (1-eHaps(t,k));
          // d3 = d1 + d2 + (1-x) * pRgivenH2(iRead); // denom
          gammaUpdate_t(k,t,0) = gammaUpdate_t(k,t,0) + d * d1;
          gammaUpdate_t(k,t,1) = gammaUpdate_t(k,t,1) + d * (d1 + d2);
        } // end of k
      } // end of SNP in read 
    } // end of read
  } // end of method 7
  //
  //
  if((model==9) | (model==10))
  {
    for(iRead=0; iRead<=nReads-1; iRead++)
    {
      // recal that below is what is used to set each element of sampleRead
      // sampleReads.push_back(Rcpp::List::create(nU,d,phiU,pRU));
      Rcpp::List readData = as<Rcpp::List>(sampleReads[iRead]);
      int J = as<int>(readData[0]); // number of SNPs on read
      int cr = as<int>(readData[1]); // central SNP in read
      arma::ivec bqU = as<arma::ivec>(readData[2]); // bq for each SNP
      arma::ivec pRU = as<arma::ivec>(readData[3]); // position of each SNP from 0 to T-1
      d3 = pRgivenH1(iRead) / (pRgivenH1(iRead) + pRgivenH2(iRead));
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
          for(k=0; k<=K-1;k++) {
              a1 = pA * eHaps_t(k,t);
              a2 = pR * (1-eHaps_t(k,t)); 
              y = a1 + a2;
              b = eMatHapOri_t(k,iRead) * d3 / y;
              d1 = a1 * b;
              d2 = a2 * b;
              d = gamma_t(k,cr) / eMatHap_t(k,iRead);
              gammaUpdate_t(k,t,0) = gammaUpdate_t(k,t,0) + d * d1;
              gammaUpdate_t(k,t,1) = gammaUpdate_t(k,t,1) + d * (d1 + d2);
          } // end of k
      } // end of SNP in read 
    } // end of read
  } // end of method 9
  //
  //
  // make jUpdate
  //
  //
  next_section="Make jUpdate";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  //
  for(t=0; t<=T-2; t++)
      for(k=0; k<=K-1; k++)
          jUpdate_t(k,t) = transMatRate_t(1, t) * alphaMat_t(k,t) * betaHat_t(k,t+1) * eMatHapSNP_t(k,t+1);
  //
  //
  //
  next_section="Done";
  prev=print_times(prev, suppressOutput, prev_section, next_section);
  prev_section=next_section;
  //
  if(whatToReturn==0)
    return(wrap(Rcpp::List::create(
      Rcpp::Named("alphaHat_t") = alphaHat_t,
      Rcpp::Named("betaHat_t") = betaHat_t,
      Rcpp::Named("gammaUpdate_t") = gammaUpdate_t,
      Rcpp::Named("eMatHapOri_t") = eMatHapOri_t,
      Rcpp::Named("gamma_t") = gamma_t,
      Rcpp::Named("jUpdate_t") = jUpdate_t,
      Rcpp::Named("eMatHap_t") = eMatHap_t)));
  if(whatToReturn==1) // update
    return(wrap(Rcpp::List::create(
      Rcpp::Named("gammaUpdate_t") = gammaUpdate_t,
      Rcpp::Named("gamma_t") = gamma_t,
      Rcpp::Named("eMatHap_t") = eMatHap_t,
      Rcpp::Named("eMatHapOri_t") = eMatHapOri_t,
      Rcpp::Named("jUpdate_t") = jUpdate_t)));
  if(whatToReturn==2) // same thing here - for simplicity 
    return(wrap(Rcpp::List::create(
      Rcpp::Named("gammaUpdate_t") = gammaUpdate_t,
      Rcpp::Named("gamma_t") = gamma_t,
      Rcpp::Named("eMatHap_t") = eMatHap_t,
      Rcpp::Named("eMatHapOri_t") = eMatHapOri_t,
      Rcpp::Named("jUpdate_t") = jUpdate_t)));
  // just stick something below to avoid warning
  return(wrap(Rcpp::List::create(
      Rcpp::Named("c") = c)));
}











//' @export
// [[Rcpp::export]]
List cpp_read_reassign(arma::ivec ord, arma::ivec qnameInteger_ord, List sampleReadsRaw, int verbose) {
  // ord is 0-based original ordering
  // qnameInteger_ord is (ordered) integer representing reads
  
  // initialize
  List sampleReads;
  int curRead = qnameInteger_ord[0];
  int iReadStart = 0;
  int nRawReads = sampleReadsRaw.size();
  arma::ivec base_bq(10000);
  arma::ivec base_pos(10000); // there shouldnt be this many SNPs
  if(verbose == 1) {
    std::cout << "curRead=" << curRead << "\n";
  }
  for (int iRead = 0; iRead < nRawReads; iRead++ ) {
    
    if(verbose == 1) {
      std::cout << "iRead=" << iRead << "\n";
      std::cout << "qnameInteger_ord[iRead + 1]=" << qnameInteger_ord[iRead + 1] << "\n";
      std::cout << "curRead==" << curRead << "\n";    
    }
    
    if (qnameInteger_ord[iRead + 1] != curRead) {
      
      int nSNPsInRead = -1;
      for(int j = iReadStart; j <= iRead; j++) {

	int r = ord[j];
	
	if(verbose == 1) {
	  std::cout << "j=" << j << "\n";
	  std::cout << "r=" << r << "\n";
	}
	
	// so say first read is 0-based 0:2
        // want to take reads from ord[0:2]
	Rcpp::List readData = as<Rcpp::List>(sampleReadsRaw[r]);
	arma::ivec bqU = as<arma::ivec>(readData[2]);
	arma::ivec pRU = as<arma::ivec>(readData[3]);
	for(int k = 0; k < int(bqU.size()); k++) {
	  // std::cout << "k=" << k << "\n";
	  nSNPsInRead++;
	  base_bq[nSNPsInRead] = bqU[k];
	  base_pos[nSNPsInRead] = pRU[k];	  
	}
      }

      arma::ivec bqL = base_bq.subvec(0, nSNPsInRead);
      arma::ivec posL = base_pos.subvec(0, nSNPsInRead);
      sampleReads.push_back(Rcpp::List::create(nSNPsInRead, 0, bqL, posL));
      iReadStart = iRead + 1;
      curRead = qnameInteger_ord[iRead + 1]; // + 1
      if(verbose == 1) {
	std::cout << "reset curRead=" << curRead << "\n";
      }
      
    }
    
  }
  
  return sampleReads;
}
