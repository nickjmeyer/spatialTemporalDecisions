#!/bin/bash

## quit on error
set -e

PROJ_ROOT=.

## runSpatial
./scripts/runSpatial.sh 2>&1 | tee >(sed 's/.*\r//'> runSpatial.log)
