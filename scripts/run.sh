#!/bin/bash

PROJ_ROOT=.

EXEC_PATH=${PROJ_ROOT}/build/main

WNS=${PROJ_ROOT}/data/wns
TOY=${PROJ_ROOT}/data/toy

DRY=--nodryRun

EDGE_TO_EDGE="--edgeToEdge ${DRY}"
SPATIAL="${SPATIAL} ${DRY}"

${EXEC_PATH}/tuneGenWNS --srcDir ${PROJ_ROOT}/data/wns ${SPATIAL} ${DRY}
${EXEC_PATH}/tuneGenWNS --srcDir ${PROJ_ROOT}/data/wns --edgeToEdge ${DRY}

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
						${SPATIAL}
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/scalefree500 \
						${SPATIAL}
${EXEC_PATH}/copyParams --srcDir ${WNS} --outDir ${TOY}/scalefree1000 \
						${SPATIAL}


########################################
## tune generative models

## spatial spread
${EXEC_PATH}/tuneGen --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/tuneGen --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/tuneGen --srcDir ${WNS} ${SPATIAL}

${EXEC_PATH}/tuneGen --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/tuneGen --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/tuneGen --srcDir ${WNS} ${SPATIAL}

${EXEC_PATH}/tuneGen --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/tuneGen --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/tuneGen --srcDir ${WNS} ${SPATIAL}

## network spread
${EXEC_PATH}/tuneGen --srcDir ${WNS} ${EDGE_TO_EDGE}
${EXEC_PATH}/tuneGen --srcDir ${WNS} ${EDGE_TO_EDGE}
${EXEC_PATH}/tuneGen --srcDir ${WNS} ${EDGE_TO_EDGE}

${EXEC_PATH}/tuneGen --srcDir ${WNS} ${EDGE_TO_EDGE}
${EXEC_PATH}/tuneGen --srcDir ${WNS} ${EDGE_TO_EDGE}
${EXEC_PATH}/tuneGen --srcDir ${WNS} ${EDGE_TO_EDGE}

${EXEC_PATH}/tuneGen --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/tuneGen --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/tuneGen --srcDir ${WNS} ${SPATIAL}



########################################
## run sims {correctly specified}

## wns
${EXEC_PATH}/runM1Mles --srcDir ${WNS} ${SPATIAL}

## spatial spread
${EXEC_PATH}/runM1Mles --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/runM1Mles --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/runM1Mles --srcDir ${WNS} ${SPATIAL}

${EXEC_PATH}/runM1Mles --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/runM1Mles --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/runM1Mles --srcDir ${WNS} ${SPATIAL}

${EXEC_PATH}/runM1Mles --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/runM1Mles --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/runM1Mles --srcDir ${WNS} ${SPATIAL}

## network spread
${EXEC_PATH}/runM1Mles --srcDir ${WNS} ${EDGE_TO_EDGE}
${EXEC_PATH}/runM1Mles --srcDir ${WNS} ${EDGE_TO_EDGE}
${EXEC_PATH}/runM1Mles --srcDir ${WNS} ${EDGE_TO_EDGE}

${EXEC_PATH}/runM1Mles --srcDir ${WNS} ${EDGE_TO_EDGE}
${EXEC_PATH}/runM1Mles --srcDir ${WNS} ${EDGE_TO_EDGE}
${EXEC_PATH}/runM1Mles --srcDir ${WNS} ${EDGE_TO_EDGE}

${EXEC_PATH}/runM1Mles --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/runM1Mles --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/runM1Mles --srcDir ${WNS} ${SPATIAL}


########################################
## run sims {incorrectly specified}

## wns
${EXEC_PATH}/runM1MlesMiss --srcDir ${WNS} ${SPATIAL}

## spatial spread
${EXEC_PATH}/runM1MlesMiss --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/runM1MlesMiss --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/runM1MlesMiss --srcDir ${WNS} ${SPATIAL}

${EXEC_PATH}/runM1MlesMiss --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/runM1MlesMiss --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/runM1MlesMiss --srcDir ${WNS} ${SPATIAL}

${EXEC_PATH}/runM1MlesMiss --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/runM1MlesMiss --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/runM1MlesMiss --srcDir ${WNS} ${SPATIAL}

## network spread
${EXEC_PATH}/runM1MlesMiss --srcDir ${WNS} ${EDGE_TO_EDGE}
${EXEC_PATH}/runM1MlesMiss --srcDir ${WNS} ${EDGE_TO_EDGE}
${EXEC_PATH}/runM1MlesMiss --srcDir ${WNS} ${EDGE_TO_EDGE}

${EXEC_PATH}/runM1MlesMiss --srcDir ${WNS} ${EDGE_TO_EDGE}
${EXEC_PATH}/runM1MlesMiss --srcDir ${WNS} ${EDGE_TO_EDGE}
${EXEC_PATH}/runM1MlesMiss --srcDir ${WNS} ${EDGE_TO_EDGE}

${EXEC_PATH}/runM1MlesMiss --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/runM1MlesMiss --srcDir ${WNS} ${SPATIAL}
${EXEC_PATH}/runM1MlesMiss --srcDir ${WNS} ${SPATIAL}
