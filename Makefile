## setup

BUILDDIR=build/

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
PROGS:=$(PROGS:%=$(BUILDDIR)%)

CPP_SRC:=$(CPP_SRC:%=src/%.cpp)
CPP_OBJ:=$(CPP_SRC:src/%.cpp=$(BUILDDIR)%.o)

LIB=$(BUILDDIR)libspatialDecisionMaking.so

## test code

CPP_SRC_TEST:=$(wildcard src/test/*.cpp)
CPP_SRC_TEST:=$(notdir $(basename $(CPP_SRC_TEST)))
CPP_SRC_TEST:=$(filter-out $(BLACKLIST),$(CPP_SRC_TEST))

PROGS_TEST:=$(CPP_SRC_TEST:%=$(BUILDDIR)test/%.bin)

CPP_OBJ_TEST:=$(CPP_SRC_TEST:%=$(BUILDDIR)test/%.o)

CPP_SRC_TEST:=$(CPP_SRC_TEST:%=src/test/%)



## options

CC=g++-4.9

CPP_FLAGS=-std=c++11 -fopenmp -g3
LD_FLAGS=-Isrc -L$(BUILDDIR) -lgsl -larmadillo -shared -fPIC

## rules

all: $(BUILDDIR) $(LIB) $(PROGS)

test: $(BUILDDIR)test $(LIB) $(PROGS_TEST)

$(BUILDDIR)%.bin: src/%.cpp $(LIB)
	$(CC) $(LD_FLAGS) $(CPP_FLAGS) -l$(LIB:$(BUILDDIR)lib%.so=%) -o $@ $^
	ln -rfs $@ $(@:%.bin=%)

$(LIB): $(CPP_OBJ)
	$(CC) $(LD_FLAGS) $(CPP_FLAGS) -o $@ $^

$(BUILDDIR)%.o: src/%.cpp $(BUILDDIR)%.d
	$(CC) $(LD_FLAGS) $(CPP_FLAGS) -c $< -o $@

$(BUILDDIR)%.d: src/%.cpp
	$(CC) $(LD_FLAGS) $(CPP_FLAGS) -MM $< -MT $(@:%.d=%.o) > $@

-include $(CPP_OBJ:%.o=%.d)

clean:
	rm -rf $(BUILDDIR)
	mkdir -p $(BUILDDIR)test
