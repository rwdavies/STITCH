#!/usr/bin/env bash

set -e

# After STITCH has been built, test the command line interface
# This isn't part of unit tests as that would fail if STITCH has never been installed before

logfile=`mktemp`
script_dir=`dirname "$0"`
cd "${script_dir}"/../
export PATH=${PATH}:`pwd`/

logfile="temp.txt"
R --slave -f ./STITCH/tests/testthat/test-cli.R
