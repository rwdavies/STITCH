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

##for i_version in $(seq 0 $((${#version_list[@]} - 1)))
for i_version in 3 3
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
    ##export R_LIBS=${LOCAL_R_LIB}
    ##export R_LIBS_USER=${LOCAL_R_LIB}
    r_libs=`R --slave -e ".libPaths()[1]" | awk '{print substr($0, 5, 1000)}' | tr -d '"'`
    rm -r -f ${r_libs}/00LOCK-STITCH
    export SEQLIB_ROOT=${STITCH_HOME}/SeqLib/
    echo Installing ${STITCH_HOME}/releases/STITCH_${version}.tar.gz    
    R CMD INSTALL ${STITCH_HOME}/releases/STITCH_${version}.tar.gz
    
    d="/data/smew1/rdavies/stitch_development/outbredmice_chr19/"
    CFW_DATA_DIR="/data/smew1/rdavies/stitch_development/truth/cfw/"
    TEST_RESULTS_DIR=${STITCH_HOME}/test-results/
    N_CORES=16
    BAMLIST=${d}/bamlist.txt
    ##d="/well/myers/rwdavies/stitch_development/outbredmice_chr19/"
    ##CFW_DATA_DIR="/well/myers/rwdavies/stitch_development/truth/cfw/" 
    ##TEST_RESULTS_DIR=/well/myers/rwdavies/stitch_development/test-results/
    ##BAMLIST=${d}/bamlist.rescomp.txt
    ##N_CORES=1
    description="whole_chr_CFW_${name}"
    OUTPUTDIR=${TEST_RESULTS_DIR}/${description}/

    SUBMIT_SCRIPT=${STITCH_HOME}/benchmark-results/${description}.sh
    rm -r -f ${SUBMIT_SCRIPT}
    if [ ${interface} == "cli" ]
    then
       export CLI_FUNCTION_BUILD=${STITCH_HOME}/releases/STITCH_${version}.tar.gz
       export CLI_VERSION=${version}
       ${STITCH_HOME}/scripts/test-cli.R
       cp ${STITCH_HOME}/STITCH.R ${LOCAL_R_LIB}/STITCH.R
       echo -e "
hostname
which time
##export R_LIBS=${LOCAL_R_LIB}
##export R_LIBS_USER=${LOCAL_R_LIB}
/usr/bin/time -v \
     		     ${LOCAL_R_LIB}/STITCH.R \
     		     --outputdir=${OUTPUTDIR} \
     		     --bamlist=${BAMLIST} \
     		     --posfile=${d}/pos.chr19.txt \
     		     --genfile=${d}/gen.chr19.txt \
     		     --chr=chr19 \
     		     --nCores=${N_CORES} \
     		     --K=4 \
     		     --nGen=100 ${extra} 2>&1 | \
     	   tee ${STITCH_HOME}/benchmark-results/${description}.txt
" > ${SUBMIT_SCRIPT}

    elif [ ${interface} == "R" ]
    then
       echo -e "	
	/usr/bin/time -v \
		      R -e \"library(STITCH); STITCH(outputdir='\"${OUTPUTDIR}/\"', bamlist='\"${BAMLIST}\"', posfile='\"${d}/pos.chr19.txt\"', genfile='\"${d}/gen.chr19.txt\"', chr='chr19', nCores=${N_CORES}, K = 4, nGen = 100, tempdir = paste0(tempdir(), '/'));\" 2>&1 | \
	    tee ${STITCH_HOME}/benchmark-results/${description}.txt
" > ${SUBMIT_SCRIPT}
    else
	echo "bad version"
	exit 1
    fi

    ## fix bug in header and then bgzip
    if [ ${version} == "1.1.1" ]
    then
   echo -e "    
gunzip ${OUTPUTDIR}/stitch.chr19.vcf.gz
sed 's/FORMAT=<ID=GT:,Number=1/FORMAT=<ID=GT,Number=1/' ${OUTPUTDIR}/stitch.chr19.vcf > ${OUTPUTDIR}/stitch.chr19.vcf.temp
mv ${OUTPUTDIR}/stitch.chr19.vcf.temp ${OUTPUTDIR}/stitch.chr19.vcf
bgzip ${OUTPUTDIR}/stitch.chr19.vcf
tabix ${OUTPUTDIR}/stitch.chr19.vcf.gz
" >> ${SUBMIT_SCRIPT}
    fi
    ## fix non bgzipped file
    if [ ${version} == "1.2.5" ]
    then
   echo -e "    
gunzip ${OUTPUTDIR}/stitch.chr19.vcf.gz
bgzip ${OUTPUTDIR}/stitch.chr19.vcf
tabix ${OUTPUTDIR}/stitch.chr19.vcf.gz
" >> ${SUBMIT_SCRIPT}
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

