BUILDDIR=./build

BLACKLIST:=bayesP obsDataStats test3 test bayesPsamplesBR test2 tuneSp

PROGS:=$(shell find ./src/ -maxdepth 1 -name "*.cpp" -exec grep -l "int main" {} \;)
PROGS:=$(notdir $(basename $(PROGS)))
PROGS:=$(filter-out $(BLACKLIST),$(PROGS))
PROGS:=$(PROGS:=.bin)

CPP_SRC:=$(wildcard src/*.cpp)
CPP_SRC:=$(filter-out $(PROGS) $(BLACKLIST),$(CPP_SRC))

CPP_OBJ:=$(CPP_SRC:./src/%.cpp=./$(BUILDDIR)/%.o)

LIB=spatialDecisionMaking.so

CC=g++

all: $(PROGS)

%.bin: $(LIB) %.o

$(LIB): $(CPP_SRC:.cpp=.o)
	$(CC) $(LDFLAGS) -o $@ $^

%.o: %.cpp %.d
	$(CC) -o $@

# %.mk:
# 	cd src && make -f $@ $(MAKECMDGOALS)

# test:
# 	cd src/test && g++ -lgtest test_rankAgent.cpp
