mkfiles=runM1Mles.mk runM1MlesMiss.mk runM1MlesWns.mk runM1MlesWnsMiss.mk

all: $(mkfiles)


%.mk:
	cd src && make -f $@ $(MAKECMDGOALS)
