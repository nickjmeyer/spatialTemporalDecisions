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
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/rand100 ${EDGE}

${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/rand500 ${EDGE}

${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/rand1000 ${EDGE}

## tune params
${EXEC_PATH}/tuneGen --srcDir ${TOY}/rand100 ${EDGE}

${EXEC_PATH}/tuneGen --srcDir ${TOY}/rand500 ${EDGE}

${EXEC_PATH}/tuneGen --srcDir ${TOY}/rand1000 ${EDGE}

## run sims
${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/rand100 ${EDGE}

${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/rand500 ${EDGE}

${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/rand1000 ${EDGE}
