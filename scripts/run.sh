#!/bin/bash

## quit on error
set -e

PROJ_ROOT=.

./scripts/runTuneLog.sh

./scripts/runSpatialLog.sh

./scripts/runEdgeToEdgeLog.sh
