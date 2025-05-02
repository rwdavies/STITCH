#!/usr/bin/env bash

set -e
script_dir=`dirname "$0"`
cd "${script_dir}"/../
./scripts/install-r-dependencies.R

export PATH=${PATH}:`pwd`/

# Install bgzip and optionally samtools into dependencies
# symlink them into the STITCH directory

# only bother if library files not presen
# Curl seems always installed?
http_stem="https://github.com/samtools/samtools/releases/download/"
samv=1.21
bcftoolsv=1.21
htslibv=1.21
mkdir -p dependencies

get_url () {
    url=${1}
    one_if_curl_installed=`which curl | wc -l`
    one_if_wget_installed=`which wget | wc -l`
    if [ ${one_if_curl_installed} == 1 ]
    then
	curl -L -O ${url}
    elif [ ${one_if_wget_installed} == 1 ]
    then
	wget ${url}
    fi
}

# This is no longer required in STITCH 1.3.0
# but leave in for now, can be useful to have easy install
force_install=${1:-nope}
zero_if_samtools_not_installed=`which samtools | wc -l`
if [ $zero_if_samtools_not_installed == 0 ] || [ "$force_install" == "samtools" ]
then
    echo install samtools
    cd dependencies
    get_url "https://github.com/samtools/samtools/releases/download/${samv}/samtools-${samv}.tar.bz2"
    bzip2 -df samtools-${samv}.tar.bz2
    tar -xvf samtools-${samv}.tar
    cd samtools-${samv}
    ./configure
    make all
    cd ../../
    ## add soft link
    dir=`pwd`
    rm -f samtools
    ln -s "${dir}/dependencies/samtools-${samv}/samtools" "${dir}/samtools"
fi


if [ "`which bcftools | wc -l`" == "0" ] || [ "$force_install" == "bcftools" ]
then
    echo install bcftools
    cd dependencies
    get_url "https://github.com/samtools/bcftools/releases/download/${bcftoolsv}/bcftools-${bcftoolsv}.tar.bz2"    
    bzip2 -df bcftools-${bcftoolsv}.tar.bz2
    tar -xvf bcftools-${bcftoolsv}.tar
    cd bcftools-${bcftoolsv}
    make all
    cd ../../
    ## add soft link
    dir=`pwd`
    rm -f bcftools
    ln -s "${dir}/dependencies/bcftools-${bcftoolsv}/bcftools" "${dir}/bcftools"
fi



if [ "`which bgzip | wc -l`" == "0" ] || [ "`which tabix | wc -l`" == "0" ] || [ "$force_install" == "bgzip" ] || [ "$force_install" == "tabix" ]
then
    echo install bgzip
    cd dependencies
    get_url "https://github.com/samtools/htslib/releases/download/${htslibv}/htslib-${htslibv}.tar.bz2"        
    bzip2 -df htslib-${htslibv}.tar.bz2
    tar -xvf htslib-${htslibv}.tar
    cd htslib-${htslibv}
    ./configure
    make all
    cd ../../
    ## add soft link
    dir=`pwd`
    rm -f bgzip
    ln -s "${dir}/dependencies/htslib-${htslibv}/bgzip" "${dir}/bgzip"
    ln -s "${dir}/dependencies/htslib-${htslibv}/tabix" "${dir}/tabix"    
fi
