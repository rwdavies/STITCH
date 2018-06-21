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
    echo Acceptance test either runs all or one
    exit 1
fi

./scripts/test-with-testthat.sh ${what_to_test}
