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



CPPFLAGS = -fopenmp -Wall #-Wl,-rpath=/usr/lib64/R/library/RInside/lib/
INCLUDE = #-I/usr/include/R/ -I/usr/lib64/R/library/Rcpp/include/
#INCLUDE += -I/usr/lib64/R/library/RInside/include/
LINKS = -larmadillo -llapack -lblas -lgsl -lgslcblas
#LINKS += -L/usr/lib64/R/lib -L/usr/lib64/R/library/RInside/lib/ -lR -lRInside
HOST = $(shell hostname)
DEBUG = -g3 -ggdb
PROD = -O3 -DNDEBUG -DBOOST_UBLAS_NDEBUG -DARMA_NO_DEBUG -DNJM_DEBUG
PROF = $(DEBUG) -pg 
COMPILE_CPP = $(CC) $(CPPFLAGS)
BINARY = trainLambda
OBJECTS = $(BINARY).o 
OBJECTS += rand.o system.o model.o modelParam.o utilities.o agent.o \
	noTrtAgent.o myopicAgent.o proximalAgent.o rankAgent.o rankAgentToy.o \
	m1SgdOptim.o m2NmOptim.o m1SimpleOptim.o runner.o dataDepth.o calcCentrality.o \
	sortMerge.o mcmc.o settings.o
DEPENDS = $(patsubst %.o, %.d, $(OBJECTS))

ifeq "$(shell hostname)" "laber-lnx4.stat.ncsu.edu"
	CPPFLAGS+= -Wl,-rpath=/usr/lib64/mpich/lib/
endif

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


-include $(DEPENDS)

$(OBJECTS):%.o: %.cpp
	$(COMPILE_CPP) $(INCLUDE) -c $*.cpp $(LINKS)
	$(COMPILE_CPP) $(INCLUDE) $*.cpp -MM -MP -MF $*.d $(LINKS)

clean:
	rm -f $(BINARY) $(OBJECTS) $(DEPENDS) $(BINARY).tar.bz2
