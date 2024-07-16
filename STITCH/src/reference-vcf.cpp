// -*- compile-command: "clang-format -i reference-vcf.cpp" -*-

#include "vcfpp/vcfpp.h"
#include <Rcpp.h>
#include <string>
#include <vector>
#include <stdint.h>

using namespace Rcpp;
using namespace vcfpp;
using namespace std;

//' @export
// [[Rcpp::export]]
List get_rhb_from_vcf(std::string vcffile, std::string region, std::string samples = "-", bool is_check = false)
{
    BcfReader br(vcffile, region, samples);
    BcfRecord var(br.header);
    int nhaps = br.nsamples * 2;
    int nsnps = 0;
    NumericVector rowsum;
    IntegerVector pos;
    vector<bool> gt;
    vector<vector<bool>> X;
    vector<string> ref, alt;
    while(br.getNextVariant(var))
    {
        var.getGenotypes(gt);
        if(is_check)
        {
            if(!var.isNoneMissing() || !var.allPhased()) continue; // skip var with missing values and non-phased
        }
        int s = 0;
        for(auto g : gt) s += g;
        rowsum.push_back(s);
        X.push_back(gt);
        ref.push_back(var.REF());
        alt.push_back(var.ALT());
        pos.push_back(var.POS());
        nsnps++;
    }
    const int B = 32;
    int nGrids = (nsnps + B - 1) / B;
    IntegerMatrix rhb_t(nGrids, nhaps);
    int d32_times_bs, imax, ihap, k;
    // check rcpp_int_contract
    for(int bs = 0; bs < nGrids; bs++)
    {
        for(ihap = 0; ihap < nhaps; ihap++)
        {
            d32_times_bs = 32 * bs;
            if(bs < (nGrids - 1))
            {
                imax = 31;
            }
            else
            {
                imax = nsnps - d32_times_bs - 1; // final one!
            }
            std::uint32_t itmp = 0;
            for(k = imax; k >= 0; k--)
            {
                itmp <<= 1;
                int j = X[d32_times_bs + k][ihap];
                itmp |= j & 0x1;
            }
            rhb_t(bs, ihap) = itmp;
        }
    }

    return List::create(Named("rhb") = rhb_t, Named("pos") = pos, Named("ref") = ref, Named("alt") = alt,
                        Named("hapRowsum") = rowsum, Named("nhaps") = nhaps);
}

//' @export
// [[Rcpp::export]]
List Rcpp_get_hap_info_from_vcf(std::string vcffile,
                                double af_cutoff,
                                std::string region,
                                std::string samples = "-",
                                bool is_check = false,
                                bool verbose = false)
{
    BcfReader br(vcffile, region, samples);
    BcfRecord var(br.header);
    int nhaps = br.nsamples * 2;
    int nsnps = 0;
    int n_common_snps = 0;
    IntegerVector L;
    LogicalVector snp_is_common;
    double local_af;
    NumericVector count_af;
    vector<bool> gt;
    vector<vector<bool>> X;
    vector<vector<int>> rare_per_hap_info(nhaps);
    // vector<string> ref, alt;
    CharacterVector ref, alt;

    const int B = 32;
    int n_skipped = 0;
    int prev_pos = -1;
    int pseudo_commons = 0;
    while(br.getNextVariant(var))
    {
        var.getGenotypes(gt);
        if(!var.isSNP() || !var.isNoneMissing() || !var.allPhased())
        {
            n_skipped += 1;
            continue;
        }
        // only keep if meets conditions
        //  - bi-allelic
        //  - snp
        //  - position increased from previous site
        if(!((var.POS() - prev_pos) > 0))
        {
            n_skipped += 1;
            continue;
        }
        int s = 0;
        for(auto g : gt) s += g;
        ref.push_back(var.REF());
        alt.push_back(var.ALT());
        L.push_back(var.POS());
        local_af = double(s) / double(nhaps);
        count_af.push_back(double(s));
        if(local_af >= af_cutoff)
        {
            snp_is_common.push_back(true);
            n_common_snps++;
            X.push_back(gt);
        }
        else
        {
            snp_is_common.push_back(false);
            // store the indices here
            for(int i = 0; i != gt.size(); ++i)
            {
                if(gt[i])
                {
                    rare_per_hap_info[i].push_back(nsnps + 1); // 1-based
                }
            }
        }
        nsnps++;
        prev_pos = var.POS();
        if(snp_is_common.size() % B == 0) pseudo_commons = 0; // reset pseudo common counts
    }

    // operate on common herre
    int nGrids = (n_common_snps + B - 1) / B;
    IntegerMatrix rhb_t(nhaps, nGrids); // X above stores efficiently this way
    int d32_times_bs, imax, ihap, k;
    // check rcpp_int_contract
    for(int bs = 0; bs < nGrids; bs++)
    {
        for(ihap = 0; ihap < nhaps; ihap++)
        {
            d32_times_bs = 32 * bs;
            if(bs < (nGrids - 1))
            {
                imax = 31;
            }
            else
            {
                imax = n_common_snps - d32_times_bs - 1; // final one!
            }
            std::uint32_t itmp = 0;
            for(k = imax; k >= 0; k--)
            {
                itmp <<= 1;
                int j = X[d32_times_bs + k][ihap];
                itmp |= j & 0x1;
            }
            rhb_t(ihap, bs) = itmp;
        }
    }

    Rcpp::DataFrame pos = Rcpp::DataFrame::create(Rcpp::Named("POS") = clone(L), Rcpp::Named("REF") = clone(ref),
                                                  Rcpp::Named("ALT") = clone(alt));

    Rcpp::NumericMatrix ref_alleleCount(nsnps, 3);
    double nhapsd = double(nhaps);
    for(int is = 0; is < nsnps; is++)
    {
        ref_alleleCount(is, 0) = count_af(is);
        ref_alleleCount(is, 1) = nhapsd;
    }
    ref_alleleCount.column(2) = ref_alleleCount.column(0) / ref_alleleCount.column(1);

    return List::create(Named("pos") = pos, Named("rhb_t") = rhb_t, Named("rare_per_hap_info") = rare_per_hap_info,
                        Named("snp_is_common") = snp_is_common, Named("ref_alleleCount") = ref_alleleCount,
                        Named("n_skipped") = n_skipped);
}
