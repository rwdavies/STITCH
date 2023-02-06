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
  while (br.getNextVariant(var)) {
    if (is_check) {
      var.getGenotypes(gt);
      if (!var.isNoneMissing() || !var.allPhased())
        continue; // skip var with missing values and non-phased
    }
    nsnps++;
  }
  br.setRegion(region); // seek back to region
  const int B = 32;
  int nGrids = (nsnps + B - 1) / B;
  IntegerMatrix rhb(nhaps, nGrids);
  std::vector<uint32_t> X(nhaps);
  std::vector<bool> ihold(nhaps);
  int i, s = 0, t = 0;
  while (br.getNextVariant(var)) {
    var.getGenotypes(gt);
    if (is_check) {
      if (!var.isNoneMissing() || !var.allPhased())
        continue; // skip var with missing values and non-phased
    }
    // update current grids
    for (i = 0; i < nhaps; i++) {
      if (t % B == 0)
        ihold[i] = gt[i];
      X[i] = (X[i] << 1) | (gt[i] != 0);
    }
    s++;
    if (s % B == 0) {
      for (i = 0; i < nhaps; i++) {
        X[i] &= ~(1 << 31);
        rhb(i, t) = (int32_t)X[i];
        if (!ihold[i])
          rhb(i, t) |= 1 << 31; // set the 32 bit to 1
      }
      t++; // update next grid
      std::fill(X.begin(), X.end(), 0);
    }
  }
  if (nGrids == t + 1) {
    // padding 0s at the end
    for (i = 0; i < nhaps; i++) {
      X[i] <<= nGrids * B - nsnps;
      X[i] &= ~(1 << 31);
      rhb(i, t) = (int32_t)X[i];
      if (!ihold[i])
        rhb(i, t) |= 1 << 31;
    }
  } else if (nGrids == t) {
    std::cerr << "no need padding\n";
  } else {
    throw std::runtime_error("something wrong!");
  }

  return rhb;
}
