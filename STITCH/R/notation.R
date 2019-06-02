## some notes on notation
## _t means "transposed " i.e. rows = K, cols = SNPs/grids
## c or m suffix denotes post-S notation for only a single haplotype reference panel
## _tc means transposed and cube e.g. eHapsCurrent_tc has entries [k, iSNP, s] where s is 1 through S on number of sets of ancestral haps
## _m means matrix i.e. priorCurrent_m has entries [k, s], sigmaCurrent_m[SNP, s]


## transMatRate_tc is for [h, snp, s] where h is from 1-2 haploid 1-3 diploid

## note - these two, in c++, are approximately the same speed
## i.e. the .slice(s). (...stuff...) is very efficient
##         for(s = 0; s < S; s++) {
##             for(iSNP = 0; iSNP < nSNPs; iSNP++) {
##                 option 1 : dosage(iSNP) += arma::sum(gamma_t.col(iSNP) % cube_eHaps_t.slice(s).col(iSNP));
## option 2:  dosage(iSNP) += arma::sum(gamma_t.col(iSNP) % eHaps_input.col(iSNP));
##             }
##         }
