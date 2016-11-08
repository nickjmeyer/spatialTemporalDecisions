import re
import os
import sys

import argparse

def readLog(category,mode):
    ## read the contents of the logs

    logFiles = {"wns": {"edge": "runEdgeToEdgeWns.log",
                        "spatial": "runSpatialWns.log"},
                "toy": {"edge": "runEdgeToEdgeToy.log",
                        "spatial": "runSpatialToy.log"} }

    logDir = "../../data/logs"

    filePath = os.path.join(logDir,logFiles[category][mode])
    with open(filePath,"r") as f:
        contents = f.read()

    return contents


def parseLog(category,mode,log):
    header = re.compile("hostName: (?P<host>[a-zA-Z0-9-]+)\s*"
                        "fileName: (?P<file>[a-zA-Z0-9/.]+)\s*"
                        "srcDir: (?P<src>[a-zA-Z0-9/.]+)\s*"
                        "datDir: (?P<dat>[a-zA-Z0-9-/.]+)\s*")

    model = re.compile("Miss")

    toy = re.compile("(?:[.a-zA-Z0-9]+/)+"
                     "(?P<network>[a-zA-Z]+)"
                     "(?P<size>[0-9]+)")

    logInfo = []

    for headerInfo in header.finditer(log):
        info = {}

        ## wns or toy
        info["category"] = category

        ## network information for toy
        if category == "toy":
            match = toy.search(headerInfo.group("src"))
            info["network"] = match.group("network")
            info["size"] = match.group("size")
        else:
            info["network"] = "NA"
            info["size"] = "NA"

        ## model specification
        if model.search(headerInfo.group("file")):
            info["model"] = 0
        else:
            info["model"] = 1


        ## simulation results
        resFile = os.path.join("../..",
                               headerInfo.group("dat"),
                               "results_.txt")
        with open(resFile,"r") as f:
            ## parse csv
            resRaw = [line.strip().split(",")
                      for line in f.readlines()]

            ## remove empty lines
            resRaw = [line for line in resRaw if len(line) == 3]

            ## strip whitespace from each token
            resRaw = [[token.strip() for token in line]
                      for line in resRaw]

            ## create dictionary for results
            resInfo = {}
            for res in resRaw:
                resInfo[res[0]] = {"mean":float(res[1]),
                                   "sd":float(res[2])}
            info["res"] = resInfo


        ## add info to log vector
        logInfo.append(info)


    return logInfo


def formatLogInfo(info):
    category = info["category"]
    network = info["network"]
    size = info["size"]
    model = info["model"]

    infoStr = ""
    infoFmt = "%s, %s, %s, %d, %s, %f, %f"
    for agent in info["res"]:
        infoStr += infoFmt % (category,network,size,model,agent,
                              info["res"][agent]["mean"],
                              info["res"][agent]["sd"])
        infoStr += "\n"
    return infoStr

def writeLogInfo(handle,logInfo):
    header = "category, network, size, model, agent, mean, sd\n"
    handle.write(header)
    for info in logInfo:
        handle.write(formatLogInfo(info))


def main(dryRun,force):
    categories = ["wns","toy"]
    modes = ["edge","spatial"]

    ## parse the log files
    logInfo = []
    for category in categories:
        for mode in modes:
            contents = readLog(category,mode)
            logInfo.extend(parseLog(category,mode,contents))

    ## format results and write to file
    if dryRun:
        writeLogInfo(sys.stdout,logInfo)
    else:
        resFile = "../../data/results/simResults.csv"
        if os.path.isfile("../../data/results/simResults.csv") and not force:
            print "Results file exists.  Use force option to overwrite the file"
        else:
            with open(resFile,"w") as f:
                writeLogInfo(f,logInfo)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dryRun",dest="dryRun",action="store_true")
    ap.add_argument("--no-dryRun",dest="dryRun",action="store_false")
    ap.set_defaults(dryRun=True)

    ap.add_argument("--force",dest="force",action="store_true")
    ap.add_argument("--no-force",dest="force",action="store_false")
    ap.set_defaults(force=False)

    args = ap.parse_args()

    main(args.dryRun,args.force)
