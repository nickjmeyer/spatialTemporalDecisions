#!/bin/bash

## quit on error
set -e

PROJ_ROOT=.

EXEC_PATH=${PROJ_ROOT}/build/src/main

WNS=${PROJ_ROOT}/data/wns
TOY=${PROJ_ROOT}/data/toy

DRY=--nodryRun

SPATIAL="--noedgeToEdge ${DRY}"
EDGE="--edgeToEdge ${DRY}"


## tune params
${EXEC_PATH}/tuneGenWNS --srcDir ${WNS} ${SPATIAL}

## run sims
${EXEC_PATH}/runM1MlesWnsMiss --srcDir ${WNS} ${SPATIAL}
