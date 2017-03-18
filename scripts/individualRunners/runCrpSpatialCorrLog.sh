#!/bin/bash

## quit on error
set -e

PROJ_ROOT=.

./scripts/runCrpSpatialCorr.sh 2>&1 | tee >(sed 's/.*\r//'> ./data/logs/runCrpSpatialCorr.log)
