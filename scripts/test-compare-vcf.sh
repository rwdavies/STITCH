#!/usr/bin/env bash

# First, run profile to generate substrate
# Second, test against prepared megamuga data

set -e

script_dir=`dirname "$0"`
cd "${script_dir}"/../
export PATH=`pwd`/:${PATH}

./scripts/profile.R

if [ ! -e megamuga_2018_06_21.RData ]
then
    curl -O http://www.well.ox.ac.uk/~rwdavies/ancillary/megamuga_2018_06_21.RData
fi

./scripts/compare_vcf_to_truth.R --test-file=./test-results/profile-one-off/stitch.chr19.vcf.gz --chr=chr19 --compare-against=megamuga --mega-save-file=megamuga_2018_06_21.RData
