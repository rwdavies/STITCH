#!/usr/bin/env bash

set -e

script_dir=`dirname "$0"`
cd "${script_dir}"/../
export PATH=`pwd`/:${PATH}

if [ "${1}" == "" ]
then
    what_to_test=unit
else
    what_to_test=unit-"${1}"
    if ! [ -e "STITCH/tests/testthat/test-${what_to_test}.R" ]
    then
	echo Unit test either runs all, or specify file test-unit-{argument}.R
	exit 1
    fi
fi

./scripts/test-with-testthat.sh ${what_to_test}
