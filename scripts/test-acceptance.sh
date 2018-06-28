#!/usr/bin/env bash

set -e

script_dir=`dirname "$0"`
cd "${script_dir}"/../
export PATH=`pwd`/:${PATH}

if [ "${1}" == "" ]
then
    what_to_test=acceptance
elif [ "${1}" == "one" ]
then
    what_to_test=acceptance-one
else
    what_to_test="${1}"
    if ! [ -e "STITCH/tests/testthat/test-acceptance-${what_to_test}.R" ]
    then
	echo Acceptance test either runs all, one, or specify file test-acceptance-{argument}.R
	exit 1
    fi
fi

./scripts/test-with-testthat.sh ${what_to_test}
