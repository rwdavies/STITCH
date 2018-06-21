#!/usr/bin/env bash

set -e

script_dir=`dirname "$0"`
cd "${script_dir}"/../
export PATH=`pwd`/:${PATH}

./scripts/test-using-testthat.sh unit
