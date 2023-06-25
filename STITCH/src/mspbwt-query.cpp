#include "mspbwt.h"
#include "timer.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

List get_quilt_rhb_in_mspbwt(const QUILT_RHB & q)
{
    Rcpp::NumericMatrix ref_alleleCount(q.nsnps, 3);
    for(int is = 0; is < q.nsnps; is++)
    {
        ref_alleleCount(is, 0) = q.ac[is];
        ref_alleleCount(is, 1) = q.nhaps;
    }
    ref_alleleCount.column(2) = ref_alleleCount.column(0) / ref_alleleCount.column(1);

    // turn Int2D into matrix/array type in R to make make_rhb_t_equality happy
    IntegerMatrix rhb_t(q.nhaps, q.nGrids);
    for(int i = 0; i < q.nhaps; i++) std::copy(q.rhb_t[i].begin(), q.rhb_t[i].end(), rhb_t(i, Rcpp::_).begin());

    Rcpp::DataFrame pos =
        Rcpp::DataFrame::create(Rcpp::Named("POS") = q.pos, Rcpp::Named("REF") = q.ref, Rcpp::Named("ALT") = q.alt);

    return List::create(Named("pos") = pos, Named("rhb_t") = rhb_t, Named("rare_per_hap_info") = q.rare_per_hap_info,
                        Named("snp_is_common") = q.snp_is_common, Named("ref_alleleCount") = ref_alleleCount,
                        Named("n_skipped") = q.n_skipped);
}

//' @export
// [[Rcpp::export]]
List quilt_mspbwt_build(const std::string & binfile,
                        const std::string & vcfpanel,
                        const std::string & samples,
                        const std::string & region,
                        int nindices,
                        int mspbwtB,
                        double maf)
{
    if(mspbwtB == 32)
    {
        MSPBWT<uint32_t> msp(nindices);
        msp.is_quilt_rhb = true;
        msp.build(vcfpanel, samples, region, maf);
        msp.save(binfile);
        return get_quilt_rhb_in_mspbwt(msp.quilt);
    }
    else if(mspbwtB == 64)
    {
        MSPBWT<uint64_t> msp(nindices);
        msp.is_quilt_rhb = true;
        msp.build(vcfpanel, samples, region, maf);
        msp.save(binfile);
        return get_quilt_rhb_in_mspbwt(msp.quilt);
    }
    else if(mspbwtB == 128)
    {
        MSPBWT<unsigned __int128> msp(nindices);
        msp.is_quilt_rhb = true;
        msp.build(vcfpanel, samples, region, maf);
        msp.save(binfile);
        return get_quilt_rhb_in_mspbwt(msp.quilt);
    }
    else
    {
        throw invalid_argument("mspbwtB must be one of 32, 64 or 128\n");
    }
}

//' @export
// [[Rcpp::export]]
void mspbwt_build(const std::string & binfile,
                  const std::string & vcfpanel,
                  const std::string & samples,
                  const std::string & region,
                  int nindices,
                  int mspbwtB,
                  double maf)
{
    if(mspbwtB == 32)
    {
        MSPBWT<uint32_t> msp(nindices);
        msp.is_quilt_rhb = false;
        msp.build(vcfpanel, samples, region, maf);
        msp.save(binfile);
    }
    else if(mspbwtB == 64)
    {
        MSPBWT<uint64_t> msp(nindices);
        msp.is_quilt_rhb = false;
        msp.build(vcfpanel, samples, region, maf);
        msp.save(binfile);
    }
    else if(mspbwtB == 128)
    {
        MSPBWT<unsigned __int128> msp(nindices);
        msp.is_quilt_rhb = false;
        msp.build(vcfpanel, samples, region, maf);
        msp.save(binfile);
    }
    else
    {
        throw invalid_argument("mspbwtB must be one of 32, 64 or 128\n");
    }
}

//' @export
// [[Rcpp::export]]
SEXP mspbwt_load(const std::string & binfile, int mspbwtB)
{
    if(mspbwtB == 32)
    {
        MSPBWT<uint32_t> * msp = new MSPBWT<uint32_t>();
        msp->load(binfile);
        Rcpp::XPtr<MSPBWT<uint32_t>> xp(msp, true);
        return (xp);
    }
    else if(mspbwtB == 64)
    {
        MSPBWT<uint64_t> * msp = new MSPBWT<uint64_t>();
        msp->load(binfile);
        Rcpp::XPtr<MSPBWT<uint64_t>> xp(msp, true);
        return (xp);
    }
    else if(mspbwtB == 128)
    {
        MSPBWT<unsigned __int128> * msp = new MSPBWT<unsigned __int128>();
        msp->load(binfile);
        Rcpp::XPtr<MSPBWT<unsigned __int128>> xp(msp, true);
        return (xp);
    }
    else
    {
        throw invalid_argument("mspbwtB must be one of 1, 32, 64 or 128\n");
    }
}

//' @export
// [[Rcpp::export]]
List mspbwt_report(SEXP xp_, const IntegerVector & z, int pbwtL, int mspbwtB)
{
    Timer tm;
    tm.clock();

    vector<int> zc = as<vector<int>>(z);
    IntMapU haplens, hapends, hapnindicies;
    if(mspbwtB == 32)
    {
        Rcpp::XPtr<MSPBWT<uint32_t>> xp(xp_);
        auto zg = xp->encodezg(zc);
        xp->report_neighourings(haplens, hapends, hapnindicies, zg, pbwtL);
    }
    else if(mspbwtB == 64)
    {
        Rcpp::XPtr<MSPBWT<uint64_t>> xp(xp_);
        auto zg = xp->encodezg(zc);
        xp->report_neighourings(haplens, hapends, hapnindicies, zg, pbwtL);
    }
    else if(mspbwtB == 128)
    {
        Rcpp::XPtr<MSPBWT<unsigned __int128>> xp(xp_);
        auto zg = xp->encodezg(zc);
        xp->report_neighourings(haplens, hapends, hapnindicies, zg, pbwtL);
    }
    else
    {
        throw invalid_argument("mspbwtB must be one of 16, 32, 64 or 128\n");
    }
    int n = haplens.size();
    vector<int> haps(n), lens(n), ends(n), nindices(n);
    n = 0;
    for(auto const & h : haplens)
    {
        haps[n] = h.first + 1; // return 1-based to R
        lens[n] = h.second;
        ends[n] = hapends[h.first];
        nindices[n] = hapnindicies[h.first];
        n++;
    }

    Rcout << "elapsed time of mspbwt insert: " << tm.abstime() << " milliseconds" << endl;
    List out = List::create(Named("haps", haps), Named("lens", lens), Named("ends", ends), Named("n", nindices));

    return out;
}
