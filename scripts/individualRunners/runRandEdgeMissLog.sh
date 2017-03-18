#!/bin/bash

## quit on error
set -e

PROJ_ROOT=.

./scripts/individualRunners/runRandEdgeMissLog.sh 2>&1 | tee >(sed 's/.*\r//'> ./data/logs/runRandEdgeMissLog.log)
