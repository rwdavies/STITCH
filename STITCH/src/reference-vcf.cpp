#include "vcfpp/vcfpp.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace vcfpp;

//' @export
// [[Rcpp::export]]
IntegerMatrix get_rhb_from_vcf(const std::string &vcffile,
                               const std::string &samples,
                               const std::string &region) {

  BcfReader br(vcffile, samples, region);
  BcfRecord var(br.header);
  int nhaps = br.nsamples * 2;
  int nsnps = 0;
  std::vector<bool> gt;
  while (br.getNextVariant(var)) {
    var.getGenotypes(gt);
    if (!var.isNoneMissing() || !var.allPhased())
      continue; // skip var with missing values and non-phased
    nsnps++;
  }
  br.setRegion(region); // seek back to region
  const int B = 32;
  int nGrids = (nsnps + B - 1) / B;
  IntegerMatrix rhb(nhaps, nGrids);
  int i, s, t = 0;
  while (br.getNextVariant(var)) {
    var.getGenotypes(gt);
    if (!var.isNoneMissing() || !var.allPhased())
      continue;
    // update current grids
    for (i = 0; i < nhaps; i++)
      rhb(i, t) = (rhb(i, t) << 1) | (gt[i] != 0);
    s++;
    if (s % B == 0) {
      t++; // update next grid
    }
  }
  if (nGrids == t + 1) {
    // padding 0s at the end
    for (i = 0; i < nhaps; i++) {
      rhb(i, t) <<= nGrids * B - nsnps;
    }
  } else if (nGrids == t) {
    std::cerr << "no need padding\n";
  } else {
    throw std::runtime_error("something wrong!");
  }

  return rhb;
}
