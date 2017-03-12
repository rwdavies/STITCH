#!/usr/bin/env bash

set -e

# Run acceptance tests on installed version of STITCH
# Technically, the "unit" tests within STITCH include both unit and acceptance tests
# These are more medium scale, real world examples
# They are useful as they identify less obvious bugs, check performance
# Ideally, a more formal check would take on the output (with seeds) to
# ensure results are invariant and/or confirm results against validation

script_dir=`dirname "$0"`
cd "${script_dir}"/../

./examples/examples.R

exit 0

