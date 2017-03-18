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


## copy params
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/crp100 ${SPATIAL}

${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/crp500 ${SPATIAL}

${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/crp1000 ${SPATIAL}

## tune params
${EXEC_PATH}/tuneGen --srcDir ${TOY}/crp100 ${SPATIAL}

${EXEC_PATH}/tuneGen --srcDir ${TOY}/crp500 ${SPATIAL}

${EXEC_PATH}/tuneGen --srcDir ${TOY}/crp1000 ${SPATIAL}

## run sims
${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/crp100 ${SPATIAL}

${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/crp500 ${SPATIAL}

${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/crp1000 ${SPATIAL}
