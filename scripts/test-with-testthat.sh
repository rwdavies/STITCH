#!/usr/bin/env bash

set -e

what_to_test=$1
if [ "${what_to_test}" != "unit" ] && [ "${what_to_test}" != "acceptance" ] && [ "${what_to_test}" != "acceptance-one" ]
then
    if ! [ -e "STITCH/tests/testthat/test-acceptance-${what_to_test}.R" ]
    then
	echo Acceptance test either runs all, one, or specify file test-acceptance-{argument}.R
	exit 1
    fi
fi

# Run unit or acceptance tests on uncompiled package using devtools::check
# This is the fastest way I'm aware of for checking the tests using compiled code but not worrying about other packaging details
# This doesn't seem to return error codes hence the crude wrapper below
# There is probably a better way to do this that maintains speed

logfile=`mktemp`
script_dir=`dirname "$0"`
cd "${script_dir}"/../
export PATH=`pwd`/:${PATH}

logfile="temp.txt"
# suppressPackageStartupMessages
# --slave
R -e 'devtools::document("STITCH"); devtools::test("STITCH", filter = "'${what_to_test}'", reporter = "summary")' 2>&1 | tee ${logfile}

# somehow this gives 0 exit code on parse failure
started_if_1=`cat ${logfile} | grep 'Testing STITCH' | wc -l`
if [ "$started_if_1" -ne "1" ]
then
    exit 1
fi

failure=`cat ${logfile} | grep ^Failed | wc -l`
warnings=`cat ${logfile} | grep ^Warnings | wc -l `
errors=`cat ${logfile} | grep Error: | wc -l `

rm ${logfile}

if [ "$failure" -gt "0" ]
then
    exit 1
fi
if [ "$warnings" -gt "0" ]
then
    exit 1
fi
if [ "$errors" -gt "0" ]
then
    exit 1
fi


exit 0
