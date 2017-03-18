#!/bin/bash

## quit on error
set -e

PROJ_ROOT=.

./scripts/individualRunners/runScalefreeEdgeMissLog.sh 2>&1 | tee >(sed 's/.*\r//'> ./data/logs/runScalefreeEdgeMissLog.log)
