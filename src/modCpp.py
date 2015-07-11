#!/usr/bin/python

import re
import sys

if len(sys.argv) != 2:
    raise ValueError("Wrong number of command line arguments.")

fname = re.sub(r"[.].*$","",sys.argv[1])
print "Base file: %s" % fname

# fname = "mcmcGravity"

########################################
## header

source = open(fname+".hpp","r").read()

## add numBurn
source = re.sub(r"(int numSamples;\n)",r"\1int numBurn;\n",source)

## add hist containers
txt = re.search(r"(std::vector<double> intcp.*;\n)",source).groups()[0]
txt = re.sub(r"([,;])",r"Burn\1",txt)
source = re.sub(r"(std::vector<double> intcp.*;\n)",r"\1"+txt,source)

source = re.sub(r"((std::vector<double> ll)(;\n))",
                r"\1\2Burn\3",
                source)


## update headers to include saveBurn
source = re.sub(r"((void sample[(][^()]*)([)];\n))",
                r"\2"+r",\nconst bool saveBurn = false\3",source)


## set par using numburn
source = re.sub(r"(setPar[(]const int i)([)].*\n)",
                r"\1,const bool fromBurn = false\2",
                source)



with open(fname+".hpp","w") as w:
    w.write(source)

########################################
## source

source = open(fname+".cpp","r").read()

## add saveBurn
source = re.sub(r"(::sample[(].*(\n*.*){1,2})([)][{])",
                r"\1,\nconst bool saveBurn\3",
                source)
source = re.sub(r"(\n.*sample[(].*)([)];\n)",
                r"\1,saveBurn\2",
                source)

## save numBurn
source = re.sub(r"(samples.numSamples = numSamples - numBurn;\n)",
                r"\1samples.numBurn = numBurn;\n",
                source)

## clear hist containers
source = re.sub(r"((samples[.].*)([.]clear[(][)].*\n))",
                r"\1\2Burn\3",
                source)
## reserve space in hist containers
source = re.sub(r"((samples[.].*)([.]reserve[(].*)"+
                r"numSamples-numBurn(.*[)].*\n))",
                r"\1\2Burn\3numBurn\4\n",
                source)

## store in container
txt = re.search(r"(.*[}]([^}]*\n){2}"+
                r"(.*// save the samples.*\n)"+
                r"((.*samples[.][^.]*[.].*\n)*))",
                source)

ins = re.sub(r"(samples[.][^.]*)([.].*\n)",
             r"\1Burn\2",
             txt.groups()[3])
ins = re.sub(r",",r",\n",ins)

source = re.sub(r"(.*[}]([^}]*\n){2}"+
                r"(.*// save the samples.*\n)"+
                r"((.*samples[.][^.]*[.].*\n)*))",
                r"if(saveBurn){\n"+
                ins+
                r"}\n"+
                txt.groups()[0],
                source)

## set par using numburn
source = re.sub(r"(setPar[(]const int i)([)][{]\n)",
                r"\1,const bool fromBurn\2",
                source)

txt = re.search(r"(setPar[(].*[)][{]\n"+
                r"([^}]*\n)*"+
                r"(((.*)([.]at.*\n))*)"+
                r"([^}]*\n)*"+
                r"([}]\n))",
                source)
if txt:
    txt = txt.groups()[0]

    par = re.search(r"(.*[.]at.*\n"+
                    r"(.*\n)*"+
                    r".*[.]at.*\n)",
                    txt)
    par = par.groups()[0]

    parBurn = re.sub(r"[.]at",
                     r"Burn.at",
                     par)

    txt = txt.replace(par,
                      "if(fromBurn){\n"+
                      parBurn+
                      "}\n"+
                      "else{\n"+
                      par+
                      "}\n")

    source = re.sub(r"(setPar[(].*[)][{]\n"+
                    r"([^}]*\n)*"+
                    r"(((.*)([.]at.*\n))*)"+
                    r"([^}]*\n)*"+
                    r"([}]\n))",
                    txt,
                    source)


with open(fname+".cpp","w") as w:
    w.write(source)
