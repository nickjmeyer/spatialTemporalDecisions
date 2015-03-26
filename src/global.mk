CC=g++
CPPFLAGS = -std=c++11 -fopenmp -Wall
INCLUDE = -I/usr/include/superlu/
LINKS = -larmadillo -llapack -lblas -lgsl -lgslcblas -lsuperlu
HOST = $(shell hostname)
DEBUG = -g3 -ggdb
PROD = -O3 -DNDEBUG -DBOOST_UBLAS_NDEBUG -DARMA_NO_DEBUG -DNJM_NO_DEBUG
PROF = $(DEBUG) -pg 
OBJECTS = $(BINARY).o 
OBJECTS += rand.o system.o utilities.o agent.o \
	noTrtAgent.o myopicAgent.o proximalAgent.o randomAgent.o \
	rankAgent.o \
	m1SpOptim.o \
	m2QOptim.o \
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
	modelRadius.o modelParamRadius.o \
	modelCave.o modelParamCave.o \
	mcmc.o mcmcRange.o mcmcCave.o mcmcRadius.o \
	mcmcGravityTimeInf.o \
	mcmcGravityTimeInfSq.o \
	mcmcGravityTimeInfSqrt.o \
	mcmcGravityTimeInfLog.o \
	mcmcGravityTimeInfExp.o \
	mcmcGravityTimeInfExpCaves.o \
	mcmcGravityTimeInfExpLCaves.o \
	mcmcGravityTimeInfExpRCaves.o \
	pardisoSymWrap.o \
	runner.o dataDepth.o calcCentrality.o \
	sortMerge.o settings.o timer.o
DEPENDS = $(patsubst %.o, %.d, $(OBJECTS))



ifeq ("${HOST}","nick-laptop")

MKLROOT = /opt/intel/composer_xe_2015/mkl
CPPFLAGS += -DNJM_USE_MKL

else ifeq ("${HOST}","nick-office")

MKLROOT = /opt/intel/composer_xe_2015/mkl
CPPFLAGS += -DNJM_USE_MKL

else ifeq ("${HOST}","laber-lnx2")

CC = g++-4.9
CPPFLAGS+= -DRANDOM_SEED__=6

else ifeq ("${HOST}","laber-lnx3")

CC = g++-4.9
CPPFLAGS += -DRANDOM_SEED__=7
MKLROOT = /opt/intel/composer_xe_2015/mkl
CPPFLAGS += -DNJM_USE_MKL

else ifeq "${HOST}" "laber-lnx4.stat.ncsu.edu"

CC = /usr/local/gcc-4.9.2/bin/g++
CPPFLAGS += -Wl,-rpath=/usr/lib64/mpich/lib/
CPPFLAGS += -DRANDOM_SEED__=8
MKLROOT = /opt/intel/composer_xe_2015/mkl
CPPFLAGS += -DNJM_USE_MKL

else ifeq "${HOST}" "opal3.stat.ncsu.edu"

CC = /usr/local/gcc-4.9.2/bin/g++
CPPFLAGS += -Wl,-rpath=/usr/lib64/mpich/lib/
CPPFLAGS += -DRANDOM_SEED__=9

endif


ifdef MKLROOT

INCLUDE +=  -m64 -I${MKLROOT}/include
LINKS += -Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64	\
-lmkl_core -lmkl_sequential -lpthread -lm

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
