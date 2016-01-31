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
	densityEst.o \
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
	rankAgent.o \
	nullOptim.o \
	allAgent.o \
	m1SpOptim.o \
	features.o \
	toyFeatures5.o \
	wnsFeatures3.o \
	param.o \
	paramIntercept.o \
	paramBeta.o \
	paramBeta2.o \
	paramGravityGDist.o \
	paramGravPowGDist.o \
	paramGDist.o \
	paramTrt.o \
	model.o \
	modelGravityGDist.o \
	model2GravityGDist.o \
	model2GPowGDist.o \
	modelGDist.o \
	mcmcGravity.o \
	mcmcGravity2.o \
	runner.o dataDepth.o calcCentrality.o \
	sortMerge.o settings.o timer.o


DEPENDS = $(patsubst %.o, %.d, $(OBJECTS))



ifeq ("${HOST}","nick-laptop")

CC = g++-4.9

# MKLROOT = /opt/intel/composer_xe_2015/mkl
# CPPFLAGS += -DNJM_USE_MKL

else ifeq ("${HOST}","nick-office")

CC = g++-4.9

# MKLROOT = /opt/intel/composer_xe_2015/mkl
# CPPFLAGS += -DNJM_USE_MKL

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
