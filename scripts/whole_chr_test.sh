#!/usr/bin/env bash

set -e

# This script runs the CFW mice for a whole chromosome, for a better check that any changes to the code have not changed the accuracy of STITCH
# This requires ~30G of local data (BAMs) and is not for general use

script_dir=`dirname "$0"`
cd "${script_dir}"/../
STITCH_HOME=`pwd`
export PATH=`pwd`/:${PATH}
mkdir -p benchmark-results


#1.2.5 1.2.9 1.3.2 1.3.4
for version in 1.3.4 1.3.5
do

    if [ 1 == 0 ]
    then
    r_libs=`R --slave -e ".libPaths()[1]" | awk '{print substr($0, 5, 1000)}' | tr -d '"'`
    rm -r -f ${r_libs}/00LOCK-STITCH
    echo ${STITCH_HOME}/releases/STITCH_${version}.tar.gz
    R CMD INSTALL ${STITCH_HOME}/releases/STITCH_${version}.tar.gz
    ${STITCH_HOME}/scripts/test-cli.R

    d="/data/smew1/rdavies/stitch_development/outbredmice_chr19/"
    /usr/bin/time -v \
		  STITCH.R \
		  --outputdir=${STITCH_HOME}/test-results/whole_chr_CFW_${version}/ \
		  --bamlist=${d}/bamlist.txt \
		  --posfile=${d}/pos.chr19.txt \
		  --genfile=${d}/gen.chr19.txt \
		  --chr=chr19 \
		  --nCores=16 \
		  --K=4 \
		  --nGen=100 2>&1 | \
	tee ${STITCH_HOME}/benchmark-results/whole_chr_CFW_${version}.txt
    fi

    ./scripts/compare_vcf_to_truth.R --vcf=${STITCH_HOME}/test-results/whole_chr_CFW_${version}/stitch.chr19.vcf.gz --chr=chr19 --compare-against=megamuga 2>&1 | \
	tee ${STITCH_HOME}/benchmark-results/whole_chr_CFW_${version}.megamuga.txt
    ./scripts/compare_vcf_to_truth.R --vcf=${STITCH_HOME}/test-results/whole_chr_CFW_${version}/stitch.chr19.vcf.gz --chr=chr19 --compare-against=affy 2>&1 | \
	tee ${STITCH_HOME}/benchmark-results/whole_chr_CFW_${version}.affy.txt

done

exit 0
