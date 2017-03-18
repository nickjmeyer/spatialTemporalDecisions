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
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/rand100 ${SPATIAL}

${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/rand500 ${SPATIAL}

${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/rand1000 ${SPATIAL}

## tune params
${EXEC_PATH}/tuneGen --srcDir ${TOY}/rand100 ${SPATIAL}

${EXEC_PATH}/tuneGen --srcDir ${TOY}/rand500 ${SPATIAL}

${EXEC_PATH}/tuneGen --srcDir ${TOY}/rand1000 ${SPATIAL}

## run sims
${EXEC_PATH}/runM1Mles --srcDir ${TOY}/rand100 ${SPATIAL}

${EXEC_PATH}/runM1Mles --srcDir ${TOY}/rand500 ${SPATIAL}

${EXEC_PATH}/runM1Mles --srcDir ${TOY}/rand1000 ${SPATIAL}
