#!/bin/bash

## quit on error
set -e

PROJ_ROOT=.

## tune
./scripts/runTune.sh 2>&1 | tee >(sed 's/.*\r//'> runTune.log)

## runSpatial
./scripts/runSpatial.sh 2>&1 | tee >(sed 's/.*\r//'> runSpatial.log)

## runEdgeToEdge
./scripts/runEdgeToEdge.sh 2>&1 | tee >(sed 's/.*\r//'> runEdgeToEdge.log)
