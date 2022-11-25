#!/usr/bin/env bash

set -e

# Run examples tests on installed version of STITCH
# STITCH has unit and acceptance tests
# These are more medium scale, real world examples
# They are useful as they identify less obvious bugs, check performance, etc
# Ideally, a more formal check would take on the output (with seeds) to
# ensure results are invariant and/or confirm results against validation

script_dir=`dirname "$0"`
cd "${script_dir}"/../

if [ "${1}" != "" ]
then
    export TEST_DIR="${1}"
fi

## remove directories, specified in the below file
rm -r -f test-results/human_tests
rm -r -f test-results/mouse_tests
./examples/examples.R

exit 0

