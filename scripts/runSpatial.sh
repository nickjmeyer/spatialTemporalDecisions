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

## wns
${EXEC_PATH}/runM1MlesWns --srcDir ${WNS} ${SPATIAL}

## spatial spread
${EXEC_PATH}/runM1Mles --srcDir ${TOY}/crp100 ${SPATIAL}
${EXEC_PATH}/runM1Mles --srcDir ${TOY}/crp500 ${SPATIAL}
${EXEC_PATH}/runM1Mles --srcDir ${TOY}/crp1000 ${SPATIAL}

${EXEC_PATH}/runM1Mles --srcDir ${TOY}/grid100 ${SPATIAL}
${EXEC_PATH}/runM1Mles --srcDir ${TOY}/grid500 ${SPATIAL}
${EXEC_PATH}/runM1Mles --srcDir ${TOY}/grid1000 ${SPATIAL}

${EXEC_PATH}/runM1Mles --srcDir ${TOY}/rand100 ${SPATIAL}
${EXEC_PATH}/runM1Mles --srcDir ${TOY}/rand500 ${SPATIAL}
${EXEC_PATH}/runM1Mles --srcDir ${TOY}/rand1000 ${SPATIAL}


########################################
## run sims {incorrectly specified}

## wns
${EXEC_PATH}/runM1MlesWnsMiss --srcDir ${WNS} ${SPATIAL}

## spatial spread
${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/crp100 ${SPATIAL}
${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/crp500 ${SPATIAL}
${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/crp1000 ${SPATIAL}

${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/grid100 ${SPATIAL}
${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/grid500 ${SPATIAL}
${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/grid1000 ${SPATIAL}

${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/rand100 ${SPATIAL}
${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/rand500 ${SPATIAL}
${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/rand1000 ${SPATIAL}
