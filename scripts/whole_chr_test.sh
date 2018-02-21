#!/usr/bin/env bash

set -e

# This script runs the CFW mice for a whole chromosome, for a better check that any changes to the code have not changed the accuracy of STITCH
# This requires ~30G of local data (BAMs) and is not for general use

run=$1
run="${run:-local}"

if [ ${run} != "local" ] && [ ${run} != "cluster" ]
then
    echo This script accepts a single argument of local or cluster
    exit 1
fi

script_dir=`dirname "$0"`
cd "${script_dir}"/../
STITCH_HOME=`pwd`
script_dir=${STITCH_HOME}/scripts/
export PATH=`pwd`/:${PATH}
mkdir -p benchmark-results

source "${script_dir}/what_to_benchmark.sh"

for i_version in $(seq 0 $((${#version_list[@]} - 1)))
do
    version=${version_list[$i_version]}
    extra=${extra_list[$i_version]}
    name=${name_list[$i_version]}
    interface=${interface_list[$i_version]}
    echo version-${version}
    echo interface-${interface}

    ## ugh, need to rebuild CLI...
    LOCAL_R_LIB=${STITCH_HOME}/test-results/whole-chr-${name}/
    mkdir -p ${LOCAL_R_LIB}
    export R_LIBS=${LOCAL_R_LIB}
    export R_LIBS_USER=${LOCAL_R_LIB}
    r_libs=`R --slave -e ".libPaths()[1]" | awk '{print substr($0, 5, 1000)}' | tr -d '"'`
    rm -r -f ${r_libs}/00LOCK-STITCH
    echo Installing ${STITCH_HOME}/releases/STITCH_${version}.tar.gz
    export SEQLIB_ROOT=${STITCH_HOME}/SeqLib/
    R CMD INSTALL ${STITCH_HOME}/releases/STITCH_${version}.tar.gz
    
    ##d="/data/smew1/rdavies/stitch_development/outbredmice_chr19/"
    ##CFW_DATA_DIR="/data/smew1/rdavies/stitch_development/truth/cfw/"    
    d="/well/myers/rwdavies/stitch_development/outbredmice_chr19/"
    CFW_DATA_DIR="/well/myers/rwdavies/stitch_development/truth/cfw/" 
    TEST_RESULTS_DIR=/well/myers/rwdavies/stitch_development/test-results/    
    description="whole_chr_CFW_${name}"
    OUTPUTDIR=${TEST_RESULTS_DIR}/${description}/

    SUBMIT_SCRIPT=${STITCH_HOME}/benchmark-results/${description}.txt
    if [ ${interface} == "cli" ]
    then
       export CLI_FUNCTION_BUILD=${STITCH_HOME}/releases/STITCH_${version}.tar.gz
       export CLI_VERSION=${version}
       ${STITCH_HOME}/scripts/test-cli.R
       cp ${STITCH_HOME}/STITCH.R ${LOCAL_R_LIB}/STITCH.R
       echo -e "
hostname
which time
export R_LIBS=${LOCAL_R_LIB}
export R_LIBS_USER=${LOCAL_R_LIB}
/usr/bin/time -v \
     		     ${LOCAL_R_LIB}/STITCH.R \
     		     --outputdir=${OUTPUTDIR} \
     		     --bamlist=${d}/bamlist.rescomp.txt \
     		     --posfile=${d}/pos.chr19.txt \
     		     --genfile=${d}/gen.chr19.txt \
     		     --chr=chr19 \
     		     --nCores=1 \
     		     --K=4 \
     		     --nGen=100 ${extra} 2>&1 | \
     	   tee ${STITCH_HOME}/benchmark-results/${description}.txt
" > ${SUBMIT_SCRIPT}

    elif [ ${interface} == "R" ]
    then
       echo -e "	
	/usr/bin/time -v \
		      R -e \"library(STITCH); STITCH(outputdir='\"${OUTPUTDIR}/\"', bamlist='\"${d}/bamlist.txt\"', posfile='\"${d}/pos.chr19.txt\"', genfile='\"${d}/gen.chr19.txt\"', chr='chr19', nCores=1, K = 4, nGen = 100, tempdir = paste0(tempdir(), '/'));\" 2>&1 | \
	    tee ${STITCH_HOME}/benchmark-results/${description}.txt
" > ${SUBMIT_SCRIPT}
    else
	echo "bad version"
	exit 1
    fi

   echo -e "    
    ## 1.1.1 require manually fixing typo in header
    ## 1.1.1, 1.2.5 require bgzipping
    ./scripts/compare_vcf_to_truth.R --vcf=${OUTPUTDIR}/stitch.chr19.vcf.gz --chr=chr19 --compare-against=megamuga --cfw-data-dir=${CFW_DATA_DIR} 2>&1 | \
	tee ${STITCH_HOME}/benchmark-results/${description}.megamuga.txt
    ./scripts/compare_vcf_to_truth.R --vcf=${OUTPUTDIR}/stitch.chr19.vcf.gz --chr=chr19 --compare-against=affy --cfw-data-dir=${CFW_DATA_DIR} 2>&1 | \
	tee ${STITCH_HOME}/benchmark-results/${description}.affy.txt
" >> ${SUBMIT_SCRIPT}
    

    if [ $run == "local" ]
    then
	bash ${SUBMIT_SCRIPT}
    else
	echo submission time!
	qsub -cwd -V -N ${name} -pe shmem 1 -q short.qc -P myers.prjc -j Y -o ${OUTPUTDIR} ${SUBMIT_SCRIPT}
    fi
    
done

exit 0
