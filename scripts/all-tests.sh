#!/usr/bin/env bash

set -e

# Run all tests while preparing release package and documentation

script_dir=`dirname "$0"`
cd "${script_dir}"/../

./scripts/install-dependencies.sh
./scripts/test-unit.sh
./scripts/test-acceptance.sh
./scripts/build-and-install.R
./scripts/test-cli.R
./scripts/test-examples.sh
