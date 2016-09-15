#!/bin/bash

PROJ_ROOT=.

EXEC_PATH=$(PROJ_ROOT)/build/main

WNS=$(PROJ_ROOT)/data/wns
TOY=$(PROJ_ROOT)/data/toy

$(EXEC_PATH)/tuneGenWNS --srcDir $(PROJ_ROOT)/data/wns --noedgeToEdge
$(EXEC_PATH)/tuneGenWNS --srcDir $(PROJ_ROOT)/data/wns --edgeToEdge

########################################
## copy params

## spatial spread
$(EXEC_PATH)/copyParams --srcDir $(WNS) --outDir $(TOY)/crp100 --noedgeToEdge
$(EXEC_PATH)/copyParams --srcDir $(WNS) --outDir $(TOY)/crp500 --noedgeToEdge
$(EXEC_PATH)/copyParams --srcDir $(WNS) --outDir $(TOY)/crp1000 --noedgeToEdge

$(EXEC_PATH)/copyParams --srcDir $(WNS) --outDir $(TOY)/grid100 --noedgeToEdge
$(EXEC_PATH)/copyParams --srcDir $(WNS) --outDir $(TOY)/grid500 --noedgeToEdge
$(EXEC_PATH)/copyParams --srcDir $(WNS) --outDir $(TOY)/grid1000 --noedgeToEdge

$(EXEC_PATH)/copyParams --srcDir $(WNS) --outDir $(TOY)/rand100 --noedgeToEdge
$(EXEC_PATH)/copyParams --srcDir $(WNS) --outDir $(TOY)/rand500 --noedgeToEdge
$(EXEC_PATH)/copyParams --srcDir $(WNS) --outDir $(TOY)/rand1000 --noedgeToEdge

## network spread
$(EXEC_PATH)/copyParams --srcDir $(WNS) --outDir $(TOY)/grid100 --edgeToEdge
$(EXEC_PATH)/copyParams --srcDir $(WNS) --outDir $(TOY)/grid500 --edgeToEdge
$(EXEC_PATH)/copyParams --srcDir $(WNS) --outDir $(TOY)/grid1000 --edgeToEdge

$(EXEC_PATH)/copyParams --srcDir $(WNS) --outDir $(TOY)/rand100 --edgeToEdge
$(EXEC_PATH)/copyParams --srcDir $(WNS) --outDir $(TOY)/rand500 --edgeToEdge
$(EXEC_PATH)/copyParams --srcDir $(WNS) --outDir $(TOY)/rand1000 --edgeToEdge

$(EXEC_PATH)/copyParams --srcDir $(WNS) --outDir $(TOY)/scalefree100 \
						--noedgeToEdge
$(EXEC_PATH)/copyParams --srcDir $(WNS) --outDir $(TOY)/scalefree500 \
						--noedgeToEdge
$(EXEC_PATH)/copyParams --srcDir $(WNS) --outDir $(TOY)/scalefree1000 \
						--noedgeToEdge


########################################
## tune generative models

## spatial spread
$(EXEC_PATH)/tuneGen --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/tuneGen --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/tuneGen --srcDir $(WNS) --noedgeToEdge

$(EXEC_PATH)/tuneGen --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/tuneGen --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/tuneGen --srcDir $(WNS) --noedgeToEdge

$(EXEC_PATH)/tuneGen --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/tuneGen --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/tuneGen --srcDir $(WNS) --noedgeToEdge

## network spread
$(EXEC_PATH)/tuneGen --srcDir $(WNS) --edgeToEdge
$(EXEC_PATH)/tuneGen --srcDir $(WNS) --edgeToEdge
$(EXEC_PATH)/tuneGen --srcDir $(WNS) --edgeToEdge

$(EXEC_PATH)/tuneGen --srcDir $(WNS) --edgeToEdge
$(EXEC_PATH)/tuneGen --srcDir $(WNS) --edgeToEdge
$(EXEC_PATH)/tuneGen --srcDir $(WNS) --edgeToEdge

$(EXEC_PATH)/tuneGen --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/tuneGen --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/tuneGen --srcDir $(WNS) --noedgeToEdge



########################################
## run sims (correctly specified)

## wns
$(EXEC_PATH)/runM1Mles --srcDir $(WNS) --noedgeToEdge

## spatial spread
$(EXEC_PATH)/runM1Mles --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/runM1Mles --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/runM1Mles --srcDir $(WNS) --noedgeToEdge

$(EXEC_PATH)/runM1Mles --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/runM1Mles --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/runM1Mles --srcDir $(WNS) --noedgeToEdge

$(EXEC_PATH)/runM1Mles --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/runM1Mles --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/runM1Mles --srcDir $(WNS) --noedgeToEdge

## network spread
$(EXEC_PATH)/runM1Mles --srcDir $(WNS) --edgeToEdge
$(EXEC_PATH)/runM1Mles --srcDir $(WNS) --edgeToEdge
$(EXEC_PATH)/runM1Mles --srcDir $(WNS) --edgeToEdge

$(EXEC_PATH)/runM1Mles --srcDir $(WNS) --edgeToEdge
$(EXEC_PATH)/runM1Mles --srcDir $(WNS) --edgeToEdge
$(EXEC_PATH)/runM1Mles --srcDir $(WNS) --edgeToEdge

$(EXEC_PATH)/runM1Mles --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/runM1Mles --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/runM1Mles --srcDir $(WNS) --noedgeToEdge


########################################
## run sims (incorrectly specified)

## wns
$(EXEC_PATH)/runM1MlesMiss --srcDir $(WNS) --noedgeToEdge

## spatial spread
$(EXEC_PATH)/runM1MlesMiss --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/runM1MlesMiss --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/runM1MlesMiss --srcDir $(WNS) --noedgeToEdge

$(EXEC_PATH)/runM1MlesMiss --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/runM1MlesMiss --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/runM1MlesMiss --srcDir $(WNS) --noedgeToEdge

$(EXEC_PATH)/runM1MlesMiss --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/runM1MlesMiss --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/runM1MlesMiss --srcDir $(WNS) --noedgeToEdge

## network spread
$(EXEC_PATH)/runM1MlesMiss --srcDir $(WNS) --edgeToEdge
$(EXEC_PATH)/runM1MlesMiss --srcDir $(WNS) --edgeToEdge
$(EXEC_PATH)/runM1MlesMiss --srcDir $(WNS) --edgeToEdge

$(EXEC_PATH)/runM1MlesMiss --srcDir $(WNS) --edgeToEdge
$(EXEC_PATH)/runM1MlesMiss --srcDir $(WNS) --edgeToEdge
$(EXEC_PATH)/runM1MlesMiss --srcDir $(WNS) --edgeToEdge

$(EXEC_PATH)/runM1MlesMiss --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/runM1MlesMiss --srcDir $(WNS) --noedgeToEdge
$(EXEC_PATH)/runM1MlesMiss --srcDir $(WNS) --noedgeToEdge
