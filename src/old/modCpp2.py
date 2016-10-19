#!/usr/bin/python

import re
import sys

if len(sys.argv) != 2:
    raise ValueError("Wrong number of command line arguments.")

# fname = "mcmcGravity"

fname = re.sub(r"[.].*$","",sys.argv[1])
print "Base file: %s" % fname



########################################
## header

source = open(fname+".hpp","r").read()

## add density estimation to header
txt = re.search(r"(([#].*\n|\n)*([#].*\n))",source).groups()[0]
source = re.sub(txt,txt+"#include \"densityEst.hpp\"\n",source)

## add function to header
source = re.sub(r"(.*void setMean[(][)][;]\n)",
                r"\1  void setMode();\n",
                source)

## get name of container
className = re.search(r"class (.*Samples)[{]\n",source).groups()[0]

## get par names
pars = re.search(r"std::vector<double> (intcp,.*);\n",source).groups()[0]
pars = pars.split(",")

with open(fname+".hpp","w") as w:
    w.write(source)



########################################
## source

source = open(fname+".cpp","r").read()

## begin source
modFn = "void "+className+"::setMode(){\n"

## beta requires extra attention
betaIn = "beta" in pars
if betaIn:
    pars.remove("beta")

## add in pars except beta
for p in pars:
    modFn += "  "+p+"Set = DensityEst("+p+").max().second;\n"


## add in beta
if betaIn:
    modFn += "\n\n"
    modFn += "  int j,k;\n"
    modFn += "  betaSet.clear();\n"
    modFn += "  for(j = 0; j < numCovar; ++j){\n"
    modFn += "    std::vector<double> betaJ;\n"
    modFn += "    for(k = 0; k < numSamples; ++k){\n"
    modFn += "      betaJ.push_back(beta.at(k*numCovar + j));\n"
    modFn += "    }\n"
    modFn += "    betaSet.push_back(DensityEst(betaJ).max().second);\n"
    modFn += "  }\n"


## close function
modFn += "}\n\n\n"

source = re.sub(r"(enum([^}]|\n)*[}];\n)",r"\1\n\n\n"+modFn,source)


with open(fname+".cpp","w") as w:
    w.write(source)
