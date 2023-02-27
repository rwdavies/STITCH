#include "vcfpp/vcfpp.h"
#include <Rcpp.h>
#include <string>
#include <vector>

using namespace Rcpp;
using namespace vcfpp;
using namespace std;

//' @export
// [[Rcpp::export]]
List get_rhb_from_vcf(std::string vcffile,
                      std::string region,
                      std::string samples = "-",
                      bool is_check = false)
{
    BcfReader br(vcffile, samples, region);
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
            if(!var.isNoneMissing() || !var.allPhased())
                continue; // skip var with missing values and non-phased
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
