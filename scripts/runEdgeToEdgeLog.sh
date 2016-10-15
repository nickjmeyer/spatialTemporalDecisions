#!/bin/bash

## quit on error
set -e

PROJ_ROOT=.

## runEdgeToEdge
./scripts/runEdgeToEdge.sh 2>&1 | tee >(sed 's/.*\r//'> runEdgeToEdge.log)
