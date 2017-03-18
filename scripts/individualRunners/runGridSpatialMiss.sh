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
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/grid100 ${SPATIAL}

${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/grid500 ${SPATIAL}

${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/grid1000 ${SPATIAL}

## tune params
${EXEC_PATH}/tuneGen --srcDir ${TOY}/grid100 ${SPATIAL}

${EXEC_PATH}/tuneGen --srcDir ${TOY}/grid500 ${SPATIAL}

${EXEC_PATH}/tuneGen --srcDir ${TOY}/grid1000 ${SPATIAL}

## run sims
${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/grid100 ${SPATIAL}

${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/grid500 ${SPATIAL}

${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/grid1000 ${SPATIAL}
