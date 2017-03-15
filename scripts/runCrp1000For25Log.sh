#!/bin/bash

## quit on error
set -e

PROJ_ROOT=.

## runSpatialToy
./scripts/runCrp1000For25.sh 2>&1 | tee >(sed 's/.*\r//'> ./data/logs/runCrp1000For25.log)
