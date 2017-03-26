#!/usr/bin/env bash

set -e
script_dir=`dirname "$0"`
cd "${script_dir}"/../
export PATH=${PATH}:`pwd`/

# Install SeqLib using extra flags. Don't both

SEQ_ROOT=SeqLib
if ! [ -f ${SEQ_ROOT}/src/libseqlib.a ] || ! [ -f ${SEQ_ROOT}/bwa/libbwa.a ] || ! [ -f ${SEQ_ROOT}/fermi-lite/libfml.a ] || ! [ -f ${SEQ_ROOT}/htslib/libhts.a ]
then
    cd ${SEQ_ROOT}
    ./configure CPPFLAGS="-fPIC" LDFLAGS="-fPIC" CFLAGS="-fPIC" CXXFLAGS="-fPIC"
    # my compiler tells me I need to do this
    # can't seem to do this with configure
    sed '/^CPPFLAGS/s/CPPFLAGS =/CPPFLAGS = -fPIC/' htslib/Makefile > htslib/Makefile.temp
    mv htslib/Makefile.temp htslib/Makefile
    make
    cd ..
fi

