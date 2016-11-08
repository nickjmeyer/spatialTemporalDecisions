#!/bin/bash

## quit on error
set -e

PROJ_ROOT=.

## runSpatialToy
./scripts/runSpatialToy.sh 2>&1 | tee >(sed 's/.*\r//'> ./data/logs/runSpatialToy.log)
