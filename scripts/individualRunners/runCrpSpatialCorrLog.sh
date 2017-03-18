#!/bin/bash

## quit on error
set -e

PROJ_ROOT=.

./scripts/individualRunners/runCrpSpatialCorr.sh 2>&1 | tee >(sed 's/.*\r//'> ./data/logs/runCrpSpatialCorr.log)
