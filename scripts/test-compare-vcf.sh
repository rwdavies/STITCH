#!/usr/bin/env bash

# First, run profile to generate substrate
# Second, test against prepared megamuga data

set -e

script_dir=`dirname "$0"`
cd "${script_dir}"/../
export PATH=`pwd`/:${PATH}

if [ ! -e megamuga_2018_06_21.RData ]
then
    #curl -k -O https://www.well.ox.ac.uk/~rwdavies/megamuga_2018_06_21.RData
    curl -k -O https://www.chg.ox.ac.uk/~rwdavies/megamuga_2018_06_21.RData
fi

mkdir -p test-data/mouse_data/
cd test-data/mouse_data/
if [ ! -e STITCH_example_2016_05_10.tgz ]
then
    #curl -k -O https://www.well.ox.ac.uk/~rwdavies/ancillary/STITCH_example_2016_05_10.tgz
    curl -k -O https://www.chg.ox.ac.uk/~rwdavies/ancillary/STITCH_example_2016_05_10.tgz
fi
tar -xzf STITCH_example_2016_05_10.tgz
cd ../../

mkdir -p test-results/profile-one-off/
./scripts/profile.R

./scripts/compare_vcf_to_truth.R --test-file=./test-results/profile-one-off/stitch.chr19.vcf.gz --chr=chr19 --compare-against=megamuga --mega-save-file=megamuga_2018_06_21.RData
