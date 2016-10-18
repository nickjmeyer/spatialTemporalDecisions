#!/bin/bash

## quit on error
set -e

PROJ_ROOT=.

EXEC_PATH=${PROJ_ROOT}/build/main

WNS=${PROJ_ROOT}/data/wns
TOY=${PROJ_ROOT}/data/toy

DRY=--nodryRun

EDGE_TO_EDGE="--edgeToEdge ${DRY}"
SPATIAL="--noedgeToEdge ${DRY}"

########################################
## run sims {correctly specified}

## wns spatial spread
${EXEC_PATH}/runM1MlesWns --srcDir ${WNS} ${SPATIAL}


########################################
## run sims {incorrectly specified}

## wns spatial spread
${EXEC_PATH}/runM1MlesWnsMiss --srcDir ${WNS} ${SPATIAL}
