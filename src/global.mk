CC=g++
CPPFLAGS = -std=c++11 -Wall -Werror -Wno-error=comment -fopenmp
INCLUDE = -I/usr/include/superlu/
LINKS = -larmadillo -llapack -lblas -lgsl -lgslcblas
HOST = $(shell hostname)
DEBUG = -g3 -ggdb
PROD = -O3 -DNDEBUG -DBOOST_UBLAS_NDEBUG -DARMA_NO_DEBUG -DNJM_NO_DEBUG
PROF = $(DEBUG) -pg 
OBJECTS = $(BINARY).o 
OBJECTS += utilities.o \
	rand.o \
	runStats.o \
	system.o \
	starts.o \
	agent.o \
	data.o \
	noTrtAgent.o \
	randomAgent.o \
	myopicAgent.o \
	proximalGDistAgent.o \
	proximalEDistAgent.o \
	rankAgent.o \
	osspAgent.o \
	m1SpOptim.o \
	m1OsspOptim.o \
	m2QOptim.o \
	features.o featuresInt.o \
	toyFeatures0.o \
	toyFeatures2.o \
	toyFeatures3.o \
	toyFeatures4.o \
	toyFeatures5.o \
	wnsFeatures0.o \
	wnsFeatures1.o \
	param.o \
	paramIntercept.o \
	paramBeta.o \
	paramGravityGDist.o \
	paramGravityEDist.o \
	paramTrend.o \
	paramTrendPow.o \
	paramTrendPowCon.o \
	paramTime.o \
	paramTimeExp.o \
	paramTimeExpCaves.o \
	paramGDist.o \
	paramGDistPow.o \
	paramGDistKern.o \
	paramRadius.o \
	paramRad.o \
	paramTrt.o \
	model.o \
	modelGravityGDist.o \
	modelGravityGDistTrend.o \
	modelGravityGDistTrendPow.o \
	modelGravityGDistTrendPowCon.o \
	modelGravityEDist.o \
	modelTimeGDist.o \
	modelTimeGDistTrend.o\
	modelTimeGDistTrendPow.o \
	modelTimeGDistTrendPowCon.o \
	modelTimeExpGDist.o \
	modelTimeExpGDistTrend.o\
	modelTimeExpGDistTrendPow.o \
	modelTimeExpGDistTrendPowCon.o \
	modelTimeExpCavesGDist.o \
	modelTimeExpCavesGDistTrend.o\
	modelTimeExpCavesGDistTrendPow.o \
	modelTimeExpCavesGDistTrendPowCon.o \
	modelTimeExpCavesEDist.o \
	modelGDist.o \
	modelGDistPow.o \
	modelGDistKern.o \
	modelRadius.o \
	modelCovar.o \
	modelRad.o \
	mcmcGravity.o \
	mcmcGravityTrend.o \
	mcmcGravityTrendPow.o \
	mcmcGravityTrendPowCon.o \
	mcmcGravityTimeInf.o \
	mcmcGravityTimeInfTrend.o \
	mcmcGravityTimeInfTrendPow.o \
	mcmcGravityTimeInfTrendPowCon.o \
	mcmcGravityTimeInfExp.o \
	mcmcGravityTimeInfExpTrend.o \
	mcmcGravityTimeInfExpTrendPow.o \
	mcmcGravityTimeInfExpTrendPowCon.o \
	mcmcGravityTimeInfExpCaves.o \
	mcmcGravityTimeInfExpCavesTrend.o \
	mcmcGravityTimeInfExpCavesTrendPow.o \
	mcmcGravityTimeInfExpCavesTrendPowCon.o \
	mcmcRad.o \
	mcmcGDist.o \
	mcmcGDistPow.o \
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
# CPPFLAGS+= -DRANDOM_SEED__=6

else ifeq ("${HOST}","laber-lnx3")

CC = g++-4.9
# CPPFLAGS += -DRANDOM_SEED__=7
MKLROOT = /opt/intel/composer_xe_2015/mkl
CPPFLAGS += -DNJM_USE_MKL

else ifeq "${HOST}" "laber-lnx4.stat.ncsu.edu"

CC = /usr/local/gcc-4.9.2/bin/g++
CPPFLAGS += -Wl,-rpath=/usr/lib64/mpich/lib/
# CPPFLAGS += -DRANDOM_SEED__=8
MKLROOT = /opt/intel/composer_xe_2015/mkl
CPPFLAGS += -DNJM_USE_MKL

else ifeq "${HOST}" "opal3.stat.ncsu.edu"

CC = /usr/local/gcc-4.9.2/bin/g++
CPPFLAGS += -Wl,-rpath=/usr/lib64/mpich/lib/
# CPPFLAGS += -DRANDOM_SEED__=9

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

