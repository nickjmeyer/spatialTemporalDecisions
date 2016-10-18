#!/bin/bash

## quit on error
set -e

PROJ_ROOT=.

## runEdgeToEdgeWns
./scripts/runEdgeToEdgeWns.sh 2>&1 | tee >(sed 's/.*\r//'> runEdgeToEdgeWns.log)
