## setup
ifdef DEBUG
BUILDDIR:=.build_debug
else
BUILDDIR:=.build_release
endif


BLACKLIST:=bayesP obsDataStats test3 test bayesPsamplesBR test2 tuneSp \
toyFeatures2Multi getCov getDist isConnected mergeClusters sample \
toyFeatures0 toyFeatures1 toyFeatures2 toyFeatures3 toyFeatures4 \
toyFeatures6 toyFeatures7 wnsFeatures0 wnsFeatures1 wnsFeatures2

## main code

PROGS:=tuneGen \
tuneGenWNS \
runM1Mles \
runM1MlesMiss \
runM1MlesWns \
runM1MlesWnsMiss

CPP_SRC:=$(wildcard src/main/*.cpp)
CPP_SRC:=$(notdir $(basename $(CPP_SRC)))
CPP_SRC:=$(filter-out $(PROGS) $(BLACKLIST),$(CPP_SRC))

PROGS:=$(PROGS:%=$(BUILDDIR)/main/%.bin)
PROG_LINKS:=$(PROGS:%.bin=%)


CPP_SRC:=$(CPP_SRC:%=src/main/%.cpp)
CPP_OBJ:=$(CPP_SRC:src/main/%.cpp=$(BUILDDIR)/main/%.o)

LIB:=$(BUILDDIR)/libspatialDecisionMaking.so


## test code

CPP_SRC_TEST:=$(wildcard src/test/*.cpp)
CPP_SRC_TEST:=$(notdir $(basename $(CPP_SRC_TEST)))
CPP_SRC_TEST:=$(filter-out $(BLACKLIST),$(CPP_SRC_TEST))

PROGS_TEST:=$(CPP_SRC_TEST:%=$(BUILDDIR)/test/%.bin)
PROG_TEST_LINKS:=$(PROGS_TEST:%.bin=%)


## options

REPO_ROOT_DIRECTORY:=$(dir $(abspath $(lastword $(MAKEFILE_LIST))))


CC=g++-5

CPP_FLAGS:= -std=c++11 -Isrc/main -fPIC -fopenmp \
-DREPO_ROOT_DIRECTORY=$(REPO_ROOT_DIRECTORY)

ifdef DEBUG
CPP_FLAGS+= -ggdb
else
CPP_FLAGS+= -O3
endif

LD_FLAGS:= -L$(BUILDDIR) -lgsl -larmadillo -lgtest -lglog -lpthread \
-lboost_system -lboost_filesystem -lgit2 -lgflags

## rules

all: | $(PROGS) $(PROG_LINKS) build


test: | $(PROGS_TEST) $(PROG_TEST_LINKS) build


.PHONY: build
build:
	ln -srfn $(BUILDDIR) $@


%/.:
	mkdir -p $(@D)


$(BUILDDIR)/main/. $(BUILDDIR)/test/.: | $(BUILDDIR)/.


.SECONDEXPANSION:


$(PROG_TEST_LINKS) $(PROG_LINKS): %: | %.bin
	ln -srf $| $@


$(BUILDDIR)/%.bin: src/%.cpp $(LIB) Makefile | $$(@D)/.
	$(CC) $(CPP_FLAGS) -o $@ $< $(LD_FLAGS) $(LIB)


$(LIB): $(CPP_OBJ) Makefile | $$(@D)/.
	$(CC) $(CPP_FLAGS) -o $@ $(CPP_OBJ) $(LD_FLAGS) -shared


$(BUILDDIR)/%.o: src/%.cpp Makefile | $$(@D)/.
	$(CC) $(CPP_FLAGS) -MMD -MP -c $< -o $@ $(LD_FLAGS)


# include dependencies
-include $(CPP_OBJ:%.o=%.d)


clean:
	rm -f build
	rm -rf $(BUILDDIR)
