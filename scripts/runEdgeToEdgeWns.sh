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

## wns network spread
${EXEC_PATH}/runM1MlesWns --srcDir ${WNS} ${EDGE_TO_EDGE}


########################################
## run sims {incorrectly specified}

## wns network spread
${EXEC_PATH}/runM1MlesWnsMiss --srcDir ${WNS} ${EDGE_TO_EDGE}
