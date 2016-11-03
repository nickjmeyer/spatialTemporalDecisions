#!/bin/bash

## quit on error
set -e

PROJ_ROOT=.

EXEC_PATH=${PROJ_ROOT}/build/src/main

WNS=${PROJ_ROOT}/data/wns
TOY=${PROJ_ROOT}/data/toy

DRY=--nodryRun

EDGE_TO_EDGE="--edgeToEdge ${DRY}"
SPATIAL="--noedgeToEdge ${DRY}"

########################################
## run sims {correctly specified}

## toy network spread
${EXEC_PATH}/runM1Mles --srcDir ${TOY}/scalefree100 ${EDGE_TO_EDGE}
${EXEC_PATH}/runM1Mles --srcDir ${TOY}/scalefree500 ${EDGE_TO_EDGE}
${EXEC_PATH}/runM1Mles --srcDir ${TOY}/scalefree1000 ${EDGE_TO_EDGE}

${EXEC_PATH}/runM1Mles --srcDir ${TOY}/grid100 ${EDGE_TO_EDGE}
${EXEC_PATH}/runM1Mles --srcDir ${TOY}/grid500 ${EDGE_TO_EDGE}
${EXEC_PATH}/runM1Mles --srcDir ${TOY}/grid1000 ${EDGE_TO_EDGE}

${EXEC_PATH}/runM1Mles --srcDir ${TOY}/rand100 ${EDGE_TO_EDGE}
${EXEC_PATH}/runM1Mles --srcDir ${TOY}/rand500 ${EDGE_TO_EDGE}
${EXEC_PATH}/runM1Mles --srcDir ${TOY}/rand1000 ${EDGE_TO_EDGE}


########################################
## run sims {incorrectly specified}

## toy network spread
${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/scalefree100 ${EDGE_TO_EDGE}
${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/scalefree500 ${EDGE_TO_EDGE}
${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/scalefree1000 ${EDGE_TO_EDGE}

${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/grid100 ${EDGE_TO_EDGE}
${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/grid500 ${EDGE_TO_EDGE}
${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/grid1000 ${EDGE_TO_EDGE}

${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/rand100 ${EDGE_TO_EDGE}
${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/rand500 ${EDGE_TO_EDGE}
${EXEC_PATH}/runM1MlesMiss --srcDir ${TOY}/rand1000 ${EDGE_TO_EDGE}
