#!/bin/bash

## quit on error
set -e

PROJ_ROOT=.

./scripts/runGridSpatialMissLog.sh 2>&1 | tee >(sed 's/.*\r//'> ./data/logs/runGridSpatialMissLog.log)
