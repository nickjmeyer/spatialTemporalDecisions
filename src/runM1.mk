CC=g++
ifeq ("$(shell hostname)","laber-lnx2")
	CC=g++-4.9
endif
ifeq ("$(shell hostname)","laber-lnx3")
	CC=g++-4.9
endif
ifeq "$(shell hostname)" "laber-lnx4.stat.ncsu.edu"
	CC=/usr/local/gcc-4.9.2/bin/g++
endif



CPPFLAGS = -std=c++11 -fopenmp -Wall
INCLUDE = 
LINKS = -larmadillo -llapack -lblas -lgsl -lgslcblas
HOST = $(shell hostname)
DEBUG = -g3 -ggdb
PROD = -O3 -DNDEBUG -DBOOST_UBLAS_NDEBUG -DARMA_NO_DEBUG -DNJM_NO_DEBUG
PROF = $(DEBUG) -pg 
BINARY = runM1
OBJECTS = $(BINARY).o 
OBJECTS += rand.o system.o utilities.o agent.o \
	noTrtAgent.o myopicAgent.o proximalAgent.o randomAgent.o \
	rankAgent.o \
	m1SpOptim.o \
	features.o featuresInt.o \
	toyFeatures2.o \
	model.o modelParam.o \
	modelGravityTimeInf.o modelParamGravityTimeInf.o \
	modelGravityTimeInfSq.o modelParamGravityTimeInfSq.o \
	modelGravityTimeInfSqrt.o modelParamGravityTimeInfSqrt.o \
	modelGravityTimeInfLog.o modelParamGravityTimeInfLog.o \
	modelGravityTimeInfExp.o modelParamGravityTimeInfExp.o \
	modelGravityTimeInfExpCaves.o modelParamGravityTimeInfExpCaves.o \
	modelGravityTimeInfExpLCaves.o modelParamGravityTimeInfExpLCaves.o \
	modelGravityTimeInfExpRCaves.o modelParamGravityTimeInfExpRCaves.o \
	modelEbola.o modelParamEbola.o \
	modelRange.o modelParamRange.o \
	modelCave.o modelParamCave.o \
	mcmc.o mcmcRange.o mcmcCave.o \
	mcmcGravityTimeInf.o \
	mcmcGravityTimeInfSq.o \
	mcmcGravityTimeInfSqrt.o \
	mcmcGravityTimeInfLog.o \
	mcmcGravityTimeInfExp.o \
	mcmcGravityTimeInfExpCaves.o \
	mcmcGravityTimeInfExpLCaves.o \
	mcmcGravityTimeInfExpRCaves.o \
	runner.o dataDepth.o calcCentrality.o \
	sortMerge.o settings.o
DEPENDS = $(patsubst %.o, %.d, $(OBJECTS))

ifeq "$(shell hostname)" "laber-lnx4.stat.ncsu.edu"
	CPPFLAGS+= -Wl,-rpath=/usr/lib64/mpich/lib/
endif

## random seeds
ifeq ("$(shell hostname)","laber-lnx2")
	CPPFLAGS+= -DRANDOM_SEED__=6
endif
ifeq ("$(shell hostname)","laber-lnx3")
	CPPFLAGS+= -DRANDOM_SEED__=7
endif
ifeq "$(shell hostname)" "laber-lnx4.stat.ncsu.edu"
	CPPFLAGS+= -DRANDOM_SEED__=8
endif



COMPILE_CPP = $(CC) $(CPPFLAGS)


all: $(BINARY)

debug: CPPFLAGS += $(DEBUG)

debug: $(BINARY)

prod: CPPFLAGS += $(PROD)

prod: $(BINARY)

prof: CPPFLAGS += $(PROF)

prof: $(BINARY)

$(BINARY): $(OBJECTS)
	$(COMPILE_CPP) $(INCLUDE) -o $@ $^ $(LINKS)
	@tar -cjf $(BINARY).tar.bz2 $(BINARY).mk \
	$$(awk -F: '{gsub(/\\/,"",$$2); gsub(/{\s}+/," ",$$2); printf $$2}' \
	$(DEPENDS) | awk '{for(i=1; i<=NF; i++) print $$i}' | sort | uniq)
	mv $(BINARY) $(BINARY).tar.bz2 ../bin/


-include $(DEPENDS)

$(OBJECTS):%.o: %.cpp
	$(COMPILE_CPP) $(INCLUDE) -c $*.cpp $(LINKS)
	$(COMPILE_CPP) $(INCLUDE) $*.cpp -MM -MP -MF $*.d $(LINKS)

clean:
	rm -f $(BINARY) $(OBJECTS) $(DEPENDS) $(BINARY).tar.bz2
