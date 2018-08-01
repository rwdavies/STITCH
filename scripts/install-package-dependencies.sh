#!/usr/bin/env bash

set -e
script_dir=`dirname "$0"`
cd "${script_dir}"/../

bin_dir=`pwd`/bin
prefix=`pwd`/install/
mkdir -p ${bin_dir}
cd ${bin_dir}

mkdir -p ${prefix}
echo install to ${prefix}
gmp_v=6.1.2
mpfr_v=4.0.1
mpc_v=1.1.0

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


cd ${bin_dir}
if ! [ -e gmp-${gmp_v}.tar ]
then
    echo download gmp
    get_url https://ftp.gnu.org/gnu/gmp/gmp-${gmp_v}.tar.bz2
    bzip2 -d gmp-${gmp_v}.tar.bz2    
fi
if ! [ -e ${prefix}/lib/libgmp.so ]
then
    echo install gmp
    rm -r -f gmp-${gmp_v}
    tar -xvf gmp-${gmp_v}.tar
    cd gmp-${gmp_v}
    ./configure --enable-static --prefix=${prefix}
    make && make check && make install
fi


cd ${bin_dir}
if ! [ -e mpfr-${mpfr_v}.tar ]
then
    echo download mpfr
    get_url https://www.mpfr.org/mpfr-${mpfr_v}/mpfr-${mpfr_v}.tar.gz
    gzip -d mpfr-${mpfr_v}.tar.gz
fi
if ! [ -e ${prefix}/lib/libmpfr.so ]
then
    echo install mpfr
    rm -r -f mpfr-${mpfr_v}
    tar -xvf mpfr-${mpfr_v}.tar
    cd mpfr-${mpfr_v}
    ./configure --enable-static --prefix=${prefix} --with-gmp=${prefix}
    make && make check && make install
fi

cd ${bin_dir}
if ! [ -e mpc-${mpc_v}.tar ]
then
    echo download mpc
    get_url https://ftp.gnu.org/gnu/mpc/mpc-${mpc_v}.tar.gz
    gzip -d mpc-${mpc_v}.tar.gz
fi
if ! [ -e ${prefix}/lib/libmpc.so ]
then
    echo echo install mpc
    cd ${bin_dir}
    rm -r -f mpc-${mpc_v}
    tar -xvf mpc-${mpc_v}.tar
    cd mpc-${mpc_v}
    ./configure --enable-static --prefix=${prefix} --with-gmp=${prefix} --with-mpfr=${prefix}
    make && make check && make install
fi

cd ${prefix}/lib/
ln -s libmpc.so.3 libmpc.so.2

echo all done!
