## setup

ifdef DEBUG
BUILDDIR=.build_debug
else
BUILDDIR=.build_release
endif

BLACKLIST:=bayesP obsDataStats test3 test bayesPsamplesBR test2 tuneSp \
toyFeatures2Multi getCov getDist isConnected mergeClusters sample \
toyFeatures0 toyFeatures1 toyFeatures2 toyFeatures3 toyFeatures4 \
toyFeatures6 toyFeatures7 wnsFeatures0 wnsFeatures1 wnsFeatures2

## make code

PROGS:=$(shell find ./src/ -maxdepth 1 -name "*.cpp" -exec grep -l "int main" {} \;)
PROGS:=$(notdir $(basename $(PROGS)))
PROGS:=$(filter-out $(BLACKLIST),$(PROGS))

CPP_SRC:=$(wildcard src/*.cpp)
CPP_SRC:=$(notdir $(basename $(CPP_SRC)))
CPP_SRC:=$(filter-out $(PROGS) $(BLACKLIST),$(CPP_SRC))

PROGS:=$(PROGS:=.bin)
PROGS:=$(PROGS:%=$(BUILDDIR)/main/%)

CPP_SRC:=$(CPP_SRC:%=src/%.cpp)
CPP_OBJ:=$(CPP_SRC:src/%.cpp=$(BUILDDIR)/main/%.o)

LIB=$(BUILDDIR)/libspatialDecisionMaking.so

## test code

CPP_SRC_TEST:=$(wildcard src/test/*.cpp)
CPP_SRC_TEST:=$(notdir $(basename $(CPP_SRC_TEST)))
CPP_SRC_TEST:=$(filter-out $(BLACKLIST),$(CPP_SRC_TEST))

PROGS_TEST:=$(CPP_SRC_TEST:%=$(BUILDDIR)/test/%.bin)

CPP_OBJ_TEST:=$(CPP_SRC_TEST:%=$(BUILDDIR)/test/%.o)

CPP_SRC_TEST:=$(CPP_SRC_TEST:%=src/test/%)

## options

CC=g++-5

ifdef DEBUG
CPP_FLAGS=-std=c++11 -ggdb
else
CPP_FLAGS=-std=c++11 -O3
endif
LD_FLAGS=-Isrc -L$(BUILDDIR) -lgsl -larmadillo -fPIC -fopenmp -lboost_filesystem -lboost_system
LD_FLAGS_TEST=$(LD_FLAGS) -lgtest -lpthread

## rules

all: | build $(BUILDDIR)/main $(LIB) $(PROGS)

test: | build $(BUILDDIR)/test $(LIB) $(PROGS_TEST)

build: $(BUILDDIR)/main $(BUILDDIR)/test
	ln -rfs $(BUILDDIR) build

$(BUILDDIR)/test:
	mkdir -p $(BUILDDIR)/test

$(BUILDDIR)/main:
	mkdir -p $(BUILDDIR)/main

$(BUILDDIR)/main/%.bin: src/%.cpp $(LIB) | build
	$(CC) $(CPP_FLAGS) -o $@ $< $(LD_FLAGS) -l$(LIB:$(BUILDDIR)/lib%.so=%)
	ln -rfs $@ $(@:%.bin=%)

$(BUILDDIR)/test/%.bin: src/test/%.cpp $(LIB) | build
	$(CC) $(CPP_FLAGS) -o $@ $< $(LD_FLAGS_TEST) -l$(LIB:$(BUILDDIR)/lib%.so=%)
	ln -rfs $@ $(@:%.bin=%)

$(LIB): $(CPP_OBJ)
	$(CC) $(CPP_FLAGS) -o $@ $^ $(LD_FLAGS) -shared

$(BUILDDIR)/main/%.o: src/%.cpp $(BUILDDIR)/main/%.d | build
	$(CC) $(CPP_FLAGS) -c $< -o $@ $(LD_FLAGS)

$(BUILDDIR)/test/%.o: src/%.cpp $(BUILDDIR)/test/%.d | build
	$(CC) $(CPP_FLAGS) -c $< -o $@ $(LD_FLAGS_TEST)

$(BUILDDIR)/main/%.d: src/%.cpp | build
	$(CC) $(CPP_FLAGS) -MM $< -MT $(@:%.d=%.o) $(LD_FLAGS) > $@

$(BUILDDIR)/test/%.d: src/%.cpp | build
	$(CC) $(CPP_FLAGS) -MM $< -MT $(@:%.d=%.o) $(LD_FLAGS_TEST) > $@

ifneq ($(MAKECMDGOALS),clean)
-include $(CPP_OBJ:%.o=%.d)
-include $(CPP_OBJ_TEST:%.o=%.d)
endif

clean:
	rm -f build
	rm -rf $(BUILDDIR)
