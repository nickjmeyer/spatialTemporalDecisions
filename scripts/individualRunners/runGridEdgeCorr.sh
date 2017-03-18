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
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/grid100 ${EDGE}

${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/grid500 ${EDGE}

${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/grid1000 ${EDGE}

## tune params
${EXEC_PATH}/tuneGen --srcDir ${TOY}/grid100 ${EDGE}

${EXEC_PATH}/tuneGen --srcDir ${TOY}/grid500 ${EDGE}

${EXEC_PATH}/tuneGen --srcDir ${TOY}/grid1000 ${EDGE}

## run sims
${EXEC_PATH}/runM1Mles --srcDir ${TOY}/grid100 ${EDGE}

${EXEC_PATH}/runM1Mles --srcDir ${TOY}/grid500 ${EDGE}

${EXEC_PATH}/runM1Mles --srcDir ${TOY}/grid1000 ${EDGE}
