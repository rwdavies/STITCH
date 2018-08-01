#!/usr/bin/env bash

set -e
script_dir=`dirname "$0"`
cd "${script_dir}"/../

bin_dir="`pwd`/bin"
prefix="`pwd`/install/"
mkdir -p "${bin_dir}"
cd "${bin_dir}"

mkdir -p "${prefix}"
echo install to "${prefix}"
xz_v=5.2.4

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


cd "${bin_dir}"
if ! [ -e xz-${xz_v}.tar ]
then
    echo download xz
    get_url https://tukaani.org/xz/xz-${xz_v}.tar.gz
    gzip -d xz-${xz_v}.tar.gz
fi
if ! [ -e ${prefix}/lib/XXX ]
then
    echo install xz
    cd "${bin_dir}"
    rm -r -f xz-${xz_v}
    tar -xvf xz-${xz_v}.tar
    cd xz-${xz_v}
    ./configure --enable-static --prefix=${prefix}
    make && make check && make install
fi
