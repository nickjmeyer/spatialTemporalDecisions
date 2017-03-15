#!/bin/bash

## quit on error
set -e

PROJ_ROOT=.

EXEC_PATH=${PROJ_ROOT}/build/src/main

WNS=${PROJ_ROOT}/data/wns
TOY=${PROJ_ROOT}/data/toy

DRY=--nodryRun

SPATIAL="--noedgeToEdge ${DRY}"

########################################
## run sims {correctly specified}

## toy spatial spread
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/crp1000 ${SPATIAL}

${EXEC_PATH}/tuneGen --srcDir ${TOY}/crp1000 ${SPATIAL}

${EXEC_PATH}/runM1Mles --srcDir ${TOY}/crp1000 ${SPATIAL}
