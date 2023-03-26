#include "mspbwt/SyllablePBWT.h"
#include "utils/timer.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//' @export
// [[Rcpp::export]]
void sypbwt_build(const std::string & binfile,
                  const std::string & vcfpanel,
                  const std::string & samples,
                  const std::string & region,
                  int nindices,
                  int mspbwtB,
                  double maf)
{
    if(mspbwtB == 64)
    {
        SyllablePBWT<unsigned long long> msp;
        msp.build(vcfpanel, samples, region);
        msp.save(binfile);
    }
    else if(mspbwtB == 128)
    {
        SyllablePBWT<unsigned __int128> msp;
        msp.build(vcfpanel, samples, region);
        msp.save(binfile);
    }
    else
    {
        throw invalid_argument("mspbwtB must be one of 64 or 128\n");
    }
}

//' @export
// [[Rcpp::export]]
SEXP sypbwt_load(const std::string & binfile, int mspbwtB)
{
    if(mspbwtB == 64)
    {
        SyllablePBWT<unsigned long long> * msp = new SyllablePBWT<unsigned long long>();
        msp->load(binfile);
        Rcpp::XPtr<SyllablePBWT<unsigned long long>> xp(msp, true);
        return (xp);
    }
    else if(mspbwtB == 128)
    {
        SyllablePBWT<unsigned __int128> * msp = new SyllablePBWT<unsigned __int128>();
        msp->load(binfile);
        Rcpp::XPtr<SyllablePBWT<unsigned __int128>> xp(msp, true);
        return (xp);
    }
    else
    {
        throw invalid_argument("mspbwtB must be one of 64 or 128\n");
    }
}

//' @export
// [[Rcpp::export]]
List sypbwt_report(SEXP xp_, const IntegerVector & z, int mspbwtL, int mspbwtB)
{
    Timer tm;
    tm.clock();

    vector<int> zc = as<vector<int>>(z);
    IntMapU haplens, hapends;
    if(mspbwtB == 64)
    {
        Rcpp::XPtr<SyllablePBWT<unsigned long long>> xp(xp_);
        auto zg = xp->encode_zg(zc);
        xp->query_zg(zg, mspbwtL, haplens, hapends);
    }
    else if(mspbwtB == 128)
    {
        Rcpp::XPtr<SyllablePBWT<unsigned __int128>> xp(xp_);
        auto zg = xp->encode_zg(zc);
        xp->query_zg(zg, mspbwtL, haplens, hapends);
    }
    else
    {
        throw invalid_argument("mspbwtB must be one of 16, 32, 64 or 128\n");
    }
    int n = haplens.size();
    vector<int> haps(n), lens(n), ends(n);
    int i = 0;
    for(auto const & h : haplens)
    {
        haps[i] = h.first + 1; // return 1-based to R
        lens[i] = h.second;
        ends[i] = hapends[h.first];
        i++;
    }

    Rcout << "elapsed time of mspbwt insert: " << tm.abstime() << " milliseconds" << endl;
    List out = List::create(Named("haps", haps), Named("lens", lens), Named("ends", ends));

    return out;
}
