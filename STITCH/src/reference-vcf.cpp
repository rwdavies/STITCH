#include "vcfpp/vcfpp.h"
#include <Rcpp.h>
#include <vector>

using namespace Rcpp;
using namespace vcfpp;

//' @export
// [[Rcpp::export]]
IntegerMatrix get_rhb_from_vcf(std::string vcffile, std::string region,
                               std::string samples = "-",
                               bool is_check = false) {

  BcfReader br(vcffile, samples, region);
  BcfRecord var(br.header);
  int nhaps = br.nsamples * 2;
  int nsnps = 0;
  std::vector<bool> gt;
  std::vector<std::vector<bool>> X;
  while (br.getNextVariant(var)) {
    var.getGenotypes(gt);
    if (is_check) {
      if (!var.isNoneMissing() || !var.allPhased())
        continue; // skip var with missing values and non-phased
    }
    X.push_back(gt);
    nsnps++;
  }
  const int B = 32;
  int nGrids = (nsnps + B - 1) / B;
  IntegerMatrix rhb_t(nGrids, nhaps);
  int d32_times_bs, imax, ihap, k;

  // check rcpp_int_contract
  for (int bs = 0; bs < nGrids; bs++) {
    for (ihap = 0; ihap < nhaps; ihap++) {
      d32_times_bs = 32 * bs;
      if (bs < (nGrids - 1)) {
        imax = 31;
      } else {
        // final one!
        imax = nsnps - d32_times_bs - 1;
      }
      std::uint32_t itmp = 0;
      for (k = imax; k >= 0; k--) {
        itmp <<= 1;
        int j = X[d32_times_bs + k][ihap];
        itmp |= j & 0x1;
      }
      rhb(bs, ihap) = itmp;
    }
  }

  return rhb_t;
}
