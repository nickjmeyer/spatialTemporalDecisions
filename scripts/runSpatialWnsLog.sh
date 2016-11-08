#!/bin/bash

## quit on error
set -e

PROJ_ROOT=.

## runSpatialWns
./scripts/runSpatialWns.sh 2>&1 | tee >(sed 's/.*\r//'> ./data/logs/runSpatialWns.log)
