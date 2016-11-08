#!/bin/bash

## quit on error
set -e

PROJ_ROOT=.

## runEdgeToEdgeToy
./scripts/runEdgeToEdgeToy.sh 2>&1 | tee >(sed 's/.*\r//'> ./data/logs/runEdgeToEdgeToy.log)
