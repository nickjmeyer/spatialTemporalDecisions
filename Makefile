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

## make code

PROGS:=tuneGen \
tuneGenWns \
runM1MlesWns \
runMwMlesMiss \
runM1MlesWns \
runMWMlesWnsMiss


CPP_SRC:=$(wildcard src/*.cpp)
CPP_SRC:=$(notdir $(basename $(CPP_SRC)))
CPP_SRC:=$(filter-out $(PROGS) $(BLACKLIST),$(CPP_SRC))


CPP_SRC_TEST:=$(wildcard src/test/*.cpp)
CPP_SRC_TEST:=$(notdir $(basename $(CPP_SRC_TEST)))
CPP_SRC_TEST:=$(filter-out $(BLACKLIST),$(CPP_SRC_TEST))


PROG_LINKS:=$(addprefix $(BUILDDIR)/main/, $(PROGS))

PROGS:=$(PROGS:=.bin)
PROGS:=$(PROGS:%=$(BUILDDIR)/main/%)

CPP_SRC:=$(CPP_SRC:%=src/%.cpp)
CPP_OBJ:=$(CPP_SRC:src/%.cpp=$(BUILDDIR)/main/%.o)

LIB:=$(BUILDDIR)/libspatialDecisionMaking.so

PROGS_TEST:=$(CPP_SRC_TEST:%=$(BUILDDIR)/test/%.bin)

PROGS_TEST_LINKS:=$(addprefix $(BUILDDIR)/test/, $(CPP_SRC_TEST))

## options

CC=g++-5

CPP_FLAGS:= -std=c++11 -Isrc -fPIC -fopenmp

ifdef DEBUG
CPP_FLAGS+= -ggdb
else
CPP_FLAGS+= -O3
endif

LD_FLAGS:= -L$(BUILDDIR) -lgsl -larmadillo
LD_FLAGS_TEST:=$(LD_FLAGS) -lgtest -lpthread

## rules

all: | $(PROGS) $(PROG_LINKS) build


test: | $(PROGS_TEST) $(PROGS_TEST_LINKS) build


.PHONY: build
build:
	ln -srf $(BUILDDIR) $@


%/.:
	mkdir -p $(@D)


$(BUILDDIR)/main/. $(BUILDDIR)/test/.: | $(BUILDDIR)/.


.SECONDEXPANSION:

$(PROGS_TEST_LINKS) $(PROG_LINKS): %: | %.bin
	ln -sr $| $@


$(BUILDDIR)/main/%.bin: src/%.cpp $(LIB) Makefile | $$(@D)/.
	$(CC) $(CPP_FLAGS) -o $@ $< $(LD_FLAGS) $(LIB)

$(BUILDDIR)/test/%.bin: src/test/%.cpp $(LIB) Makefile | $$(@D)/.
	$(CC) $(CPP_FLAGS) -o $@ $< $(LD_FLAGS_TEST) $(LIB)


$(LIB): $(CPP_OBJ) Makefile | $$(@D)/.
	$(CC) $(CPP_FLAGS) -o $@ $(CPP_OBJ) $(LD_FLAGS) -shared


$(BUILDDIR)/main/%.o: src/%.cpp Makefile | $$(@D)/.
	$(CC) $(CPP_FLAGS) -MMD -MP -c $< -o $@ $(LD_FLAGS)

$(BUILDDIR)/test/%.o: src/test/%.cpp Makefile | $$(@D)/.
	$(CC) $(CPP_FLAGS) -MMD -MP -c $< -o $@ $(LD_FLAGS_TEST)


# include dependencies
-include $(CPP_OBJ:%.o=%.d)


clean:
	rm -f build
	rm -rf $(BUILDDIR)
