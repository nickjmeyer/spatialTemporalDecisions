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
## tune wns

${EXEC_PATH}/tuneGenWNS --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/tuneGenWNS --srcDir ${WNS} ${EDGE_TO_EDGE}

########################################
## copy params

## spatial spread
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/crp100 ${SPATIAL}
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/crp500 ${SPATIAL}
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/crp1000 ${SPATIAL}

${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/grid100 ${SPATIAL}
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/grid500 ${SPATIAL}
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/grid1000 ${SPATIAL}

${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/rand100 ${SPATIAL}
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/rand500 ${SPATIAL}
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/rand1000 ${SPATIAL}

## network spread
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/grid100 ${EDGE_TO_EDGE}
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/grid500 ${EDGE_TO_EDGE}
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/grid1000 ${EDGE_TO_EDGE}

${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/rand100 ${EDGE_TO_EDGE}
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/rand500 ${EDGE_TO_EDGE}
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/rand1000 ${EDGE_TO_EDGE}

${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/scalefree100 \
						${EDGE_TO_EDGE}
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/scalefree500 \
						${EDGE_TO_EDGE}
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/scalefree1000 \
						${EDGE_TO_EDGE}


########################################
## tune generative models

## wns
## these are tuned above before copying the parameters

## spatial spread
${EXEC_PATH}/tuneGen --srcDir ${TOY}/crp100 ${SPATIAL}
${EXEC_PATH}/tuneGen --srcDir ${TOY}/crp500 ${SPATIAL}
${EXEC_PATH}/tuneGen --srcDir ${TOY}/crp1000 ${SPATIAL}

${EXEC_PATH}/tuneGen --srcDir ${TOY}/grid100 ${SPATIAL}
${EXEC_PATH}/tuneGen --srcDir ${TOY}/grid500 ${SPATIAL}
${EXEC_PATH}/tuneGen --srcDir ${TOY}/grid1000 ${SPATIAL}

${EXEC_PATH}/tuneGen --srcDir ${TOY}/rand100 ${SPATIAL}
${EXEC_PATH}/tuneGen --srcDir ${TOY}/rand500 ${SPATIAL}
${EXEC_PATH}/tuneGen --srcDir ${TOY}/rand1000 ${SPATIAL}

## network spread
${EXEC_PATH}/tuneGen --srcDir ${TOY}/scalefree100 ${EDGE_TO_EDGE}
${EXEC_PATH}/tuneGen --srcDir ${TOY}/scalefree500 ${EDGE_TO_EDGE}
${EXEC_PATH}/tuneGen --srcDir ${TOY}/scalefree1000 ${EDGE_TO_EDGE}

${EXEC_PATH}/tuneGen --srcDir ${TOY}/grid100 ${EDGE_TO_EDGE}
${EXEC_PATH}/tuneGen --srcDir ${TOY}/grid500 ${EDGE_TO_EDGE}
${EXEC_PATH}/tuneGen --srcDir ${TOY}/grid1000 ${EDGE_TO_EDGE}

${EXEC_PATH}/tuneGen --srcDir ${TOY}/rand100 ${EDGE_TO_EDGE}
${EXEC_PATH}/tuneGen --srcDir ${TOY}/rand500 ${EDGE_TO_EDGE}
${EXEC_PATH}/tuneGen --srcDir ${TOY}/rand1000 ${EDGE_TO_EDGE}
