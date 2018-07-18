#!/usr/bin/env bash

set -e

run=$1
run="${run:-local}"

if [ ${run} != "local" ] && [ ${run} != "cluster" ]
then
    echo This script accepts a single argument of local or cluster
    exit 1
fi

# Profile one version against another / multiple other versions, for overall time / RAM, and for component breakdown using Rprof

script_dir=`dirname "$0"`
cd "${script_dir}"/../
STITCH_HOME=`pwd`
export PATH=`pwd`/:${PATH}
mkdir -p benchmark-results

# argh, cannot access internet on some clusters. get local version of proftools
PROFTOOLS_TARBALL=proftools_0.99-2.tar.gz
if [ ! -f ${PROFTOOLS_TARBALL} ]
then
    wget https://cran.r-project.org/src/contrib/${PROFTOOLS_TARBALL}
fi

source "${script_dir}/what_to_benchmark.sh"

for i_version in $(seq 0 $((${#version_list[@]} - 1)))
do
    version=${version_list[$i_version]}
    extra=${extra_list[$i_version]}
    option=${options_list[$i_version]}    
    name=${name_list[$i_version]}
    extension=${extension_list[$i_version]}
    echo version=${version}
    echo extra=${extra}
    echo name=${name}
    echo option=${option}
    echo extension=${extension}
    cd ${STITCH_HOME}
    # note - not sure I need this anymore if installing locally
    # I think below should work on other machines
    # Basically want the user location for libraries
    r_libs=`R --slave -e ".libPaths()[1]" | awk '{print substr($0, 5, 1000)}' | tr -d '"'`
    rm -r -f ${r_libs}/00LOCK-STITCH
    # install copy here. cannot install properly on cluster
    LOCAL_R_LIB=${STITCH_HOME}/test-results/profile-${name}/
    mkdir -p ${LOCAL_R_LIB}
    ## if commented out, just use local system install\
    ##export R_LIBS=${LOCAL_R_LIB}
    ##export R_LIBS_USER=${LOCAL_R_LIB}
    ## R --slave -e ".libPaths()"
    echo Installing ${STITCH_HOME}/releases/STITCH_${version}.tar.gz
    export SEQLIB_ROOT=${STITCH_HOME}/SeqLib/
    R CMD INSTALL ${STITCH_HOME}/${PROFTOOLS_TARBALL}
    R CMD INSTALL ${STITCH_HOME}/releases/STITCH_${version}.tar.gz
    # try for local install
    cd ${STITCH_HOME}/test-data/mouse_data/
    for use in CRAMS BAMS
    do
	for K in 4 20
	do
	    ## do local install
	    description="benchmark_${use}_${K}_${name}"
	    export OUTPUTDIR=${STITCH_HOME}/test-results/${description}/
	    mkdir -p ${OUTPUTDIR}
	    echo -e "
            ##export R_LIBS=${LOCAL_R_LIB}
            ##export R_LIBS_USER=${LOCAL_R_LIB}
	    export USE=${use}
	    export K=${K}
	    export OPTION=${option}
	    export N_CORES=1
	    export TITLE=${use}_${version}
            export OUTPUTDIR=${OUTPUTDIR}
	    export OUTPUT_PLOT=${STITCH_HOME}/benchmark-results/${description}.pdf
            ##ulimit -n 4096 ## if no internet, something about keeping lots of reference file handles available. something I need to fix. see seqlib issue
	    /usr/bin/time -v ${STITCH_HOME}/scripts/profile.R 2>&1 | tee ${STITCH_HOME}/benchmark-results/${description}.txt
" > ${OUTPUTDIR}/script.sh
	    if [ $run == "local" ]
	    then
		bash ${OUTPUTDIR}/script.sh &
	    else
		qsub -cwd -V -N ${description} -pe shmem 1 -q short.qc -P myers.prjc -j Y -o ${OUTPUTDIR} ${OUTPUTDIR}/script.sh
	    fi
	    ## 
	    ##${STITCH_HOME}/scripts/compare_vcf_to_truth.R \
	    ##    --vcf=${OUTPUTDIR}/stitch.chr19.vcf.gz \
	    ##	--chr=chr19 \
	    ##	--compare-against=megamuga \
	    ##	--mega-save-file=${STITCH_HOME}/megamuga19.RData \
	    ##	2>&1 | \
	    ##	tee ${OUTPUTDIR}/${description}.megamuga.txt
	done
    done
    if [ $run == "local" ]
    then
	wait
    fi
done

exit 0
