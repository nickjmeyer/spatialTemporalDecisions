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
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/scalefree100 ${EDGE}

${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/scalefree500 ${EDGE}

${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/scalefree1000 ${EDGE}

## tune params
${EXEC_PATH}/tuneGen --srcDir ${TOY}/scalefree100 ${EDGE}

${EXEC_PATH}/tuneGen --srcDir ${TOY}/scalefree500 ${EDGE}

${EXEC_PATH}/tuneGen --srcDir ${TOY}/scalefree1000 ${EDGE}

## run sims
${EXEC_PATH}/runM1Mles --srcDir ${TOY}/scalefree100 ${EDGE}

${EXEC_PATH}/runM1Mles --srcDir ${TOY}/scalefree500 ${EDGE}

${EXEC_PATH}/runM1Mles --srcDir ${TOY}/scalefree1000 ${EDGE}
