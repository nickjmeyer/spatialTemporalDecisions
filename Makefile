mkfiles=runM1Mles.mk runM1MlesMiss.mk runM1MlesWns.mk runM1MlesWnsMiss.mk

all: $(mkfiles)


%.mk:
	cd src && make -f $@ $(MAKECMDGOALS)

test:
	cd src/test && g++ -lgtest test_rankAgent.cpp
