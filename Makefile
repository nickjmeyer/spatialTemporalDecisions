BUILDDIR=build/

BLACKLIST:=bayesP obsDataStats test3 test bayesPsamplesBR test2 tuneSp \
toyFeatures2Multi getCov getDist isConnected mergeClusters sample \
toyFeatures0 toyFeatures1 toyFeatures2 toyFeatures3 toyFeatures4 \
toyFeatures6 toyFeatures7 wnsFeatures0 wnsFeatures1 wnsFeatures2

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

CC=g++

CPP_FLAGS=-std=c++11 -fopenmp
LD_FLAGS=-Isrc -L$(BUILDDIR) -lgsl -larmadillo -shared -fPIC

all: $(LIB) $(PROGS)

test: $(LIB) $(CPP_SRC_TEST)

$(BUILDDIR)%.bin: src/%.cpp $(LIB)
	$(CC) $(LD_FLAGS) $(CPP_FLAGS) -l$(LIB:$(BUILDDIR)lib%.so=%) -o $@ $^
	ln -rs $@ $(@:%.bin=%)

$(LIB): $(CPP_OBJ)
	$(CC) $(LD_FLAGS) $(CPP_FLAGS) -o $@ $^

$(BUILDDIR)%.o: src/%.cpp
	$(CC) $(LD_FLAGS) $(CPP_FLAGS) -c $< -o $@

$(BUILDDIR)%.d: src/%.cpp
	$(CC) $(LD_FLAGS) $(CPP_FLAGS) -MM $< -MT $(@:%.d=%.o) > $@

-include $(CPP_OBJ:%.o=%.d)
