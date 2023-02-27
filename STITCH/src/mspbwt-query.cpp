#include "mspbwt/mspbwt.h"
#include "utils/timer.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//' @export
// [[Rcpp::export]]
void mspbwt_build(const std::string& binfile, const std::string& vcfpanel,
                  const std::string& samples, const std::string& region, int nindices,
                  int mspbwtB, double maf)
{
    if (mspbwtB == 16)
    {
        msPBWT<uint16_t> msp(nindices);
        msp.build(vcfpanel, samples, region, maf);
        msp.save(binfile);
    }
    else if (mspbwtB == 32)
    {
        msPBWT<uint32_t> msp(nindices);
        msp.build(vcfpanel, samples, region, maf);
        msp.save(binfile);
    }
    else if (mspbwtB == 64)
    {
        msPBWT<uint64_t> msp(nindices);
        msp.build(vcfpanel, samples, region, maf);
        msp.save(binfile);
    }
    else if (mspbwtB == 128)
    {
        msPBWT<unsigned __int128> msp(nindices);
        msp.build(vcfpanel, samples, region, maf);
        msp.save(binfile);
    }
    else
    {
        throw invalid_argument("mspbwtB must be one of 16, 32, 64 or 128\n");
    }
}

//' @export
// [[Rcpp::export]]
SEXP mspbwt_load(const std::string& binfile, int mspbwtB)
{
    if (mspbwtB == 16)
    {
        msPBWT<uint16_t>* msp = new msPBWT<uint16_t>();
        msp->load(binfile);
        Rcpp::XPtr<msPBWT<uint16_t>> xp(msp, true);
        return (xp);
    }
    else if (mspbwtB == 32)
    {
        msPBWT<uint32_t>* msp = new msPBWT<uint32_t>();
        msp->load(binfile);
        Rcpp::XPtr<msPBWT<uint32_t>> xp(msp, true);
        return (xp);
    }
    else if (mspbwtB == 64)
    {
        msPBWT<uint64_t>* msp = new msPBWT<uint64_t>();
        msp->load(binfile);
        Rcpp::XPtr<msPBWT<uint64_t>> xp(msp, true);
        return (xp);
    }
    else if (mspbwtB == 128)
    {
        msPBWT<unsigned __int128>* msp = new msPBWT<unsigned __int128>();
        msp->load(binfile);
        Rcpp::XPtr<msPBWT<unsigned __int128>> xp(msp, true);
        return (xp);
    }
    else
    {
        throw invalid_argument("mspbwtB must be one of 16, 32, 64 or 128\n");
    }
}


//' @export
// [[Rcpp::export]]
List mspbwt_report(SEXP xp_, const IntegerVector& z, int pbwtL, int mspbwtB)
{
    Timer tm;
    tm.clock();

    vector<int> zc = as<vector<int>>(z);
    IntMapU haplens, hapends, hapnindicies;
    if (mspbwtB == 16)
    {
        Rcpp::XPtr<msPBWT<uint16_t>> xp(xp_);
        auto zg = xp->encodezg(zc);
        xp->report_setmaximal(haplens, hapends, hapnindicies, zg);
    }
    else if (mspbwtB == 32)
    {
        Rcpp::XPtr<msPBWT<uint32_t>> xp(xp_);
        auto zg = xp->encodezg(zc);
        xp->report_setmaximal(haplens, hapends, hapnindicies, zg);
    }
    else if (mspbwtB == 64)
    {
        Rcpp::XPtr<msPBWT<uint64_t>> xp(xp_);
        auto zg = xp->encodezg(zc);
        xp->report_setmaximal(haplens, hapends, hapnindicies, zg);
    }
    else if (mspbwtB == 128)
    {
        Rcpp::XPtr<msPBWT<unsigned __int128>> xp(xp_);
        auto zg = xp->encodezg(zc);
        xp->report_setmaximal(haplens, hapends, hapnindicies, zg);
    }
    else
    {
        throw invalid_argument("mspbwtB must be one of 16, 32, 64 or 128\n");
    }
    int n = haplens.size();
    vector<int> haps(n), lens(n), ends(n), nindices(n);
    n = 0;
    for (auto const& h : haplens)
    {
        haps[n] = h.first + 1; // return 1-based to R
        lens[n] = h.second;
        ends[n] = hapends[h.first];
        nindices[n] = hapnindicies[h.first];
        n++;
    }

    Rcout << "elapsed time of mspbwt insert: " << tm.abstime() << " milliseconds" << endl;
    List out = List::create(Named("haps", haps), Named("lens", lens), Named("ends", ends),
                            Named("n", nindices));

    return out;
}
