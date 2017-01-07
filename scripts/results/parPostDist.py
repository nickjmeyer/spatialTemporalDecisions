"""Get the posterior distributions"""

import numpy as np
import os

def getPostDist(samplesFile, order, fns = {}):
    with open(samplesFile,"r") as f:
       lines = [map(float,line.strip().split()) for line in f.readlines()
                 if len(line.strip().split()) > 0]

       samples = np.array(lines)

    assert len(set(order)) == len(order)
    assert min(order) == 0
    assert max(order) == (samples.shape[1] - 1)

    postDist = []
    for i, index in enumerate(order):
        if i in fns:
            s = np.asarray(map(fns[i], samples[:,index]))
        else:
            s = samples[:,index]

        postMean = np.mean(s)
        postL95 = np.percentile(s, 2.5)
        postU95 = np.percentile(s, 97.5)

        postDist.append((postMean, postL95, postU95))

    return postDist

def getSimVals(directory,parFiles):
    simVals = []
    for pf in parFiles:
        with open(os.path.join(directory,pf)) as f:
            simVals.append(float(f.readlines()[0].strip()))

    return simVals


def main():
    ## edge
    edgeSamplesFile = ("../../data/wns/2016-11-03-16-24-26/"
                   "sampStats_edgeToEdge2_param_edge_.txt")
    edgeSamplesOrder = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    edgeDir = ("../../data/wns/Param2EdgeToEdge")
    edgeParFiles = ("intcp.txt","beta0.txt","beta1.txt","beta2.txt",
                    "beta3.txt","beta4.txt","beta5.txt","beta6.txt",
                    "beta7.txt","trtAct.txt","trtPre.txt")
    edgeParNames = ("$\\theta_0$", "$\\theta_{1,0}$", "$\\theta_{1,1}$",
                    "$\\theta_{1,2}$", "$\\theta_{1,3}$", "$\\theta_{2,0}$",
                    "$\\theta_{2,1}$", "$\\theta_{2,2}$", "$\\theta_{2,3}$",
                    "$\\theta_3$", "$\\theta_4$")

    ## spatial
    spatialSamplesFile = ("../../data/wns/2016-11-03-16-24-26/"
                   "sampStats_gravity2_param_spatial_.txt")
    spatialSamplesOrder = [0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 9, 10]
    spatialDir = ("../../data/wns/Param2GravityEDist")
    spatialParFiles = ("intcp.txt","beta0.txt","beta1.txt","beta2.txt",
                       "beta3.txt","beta4.txt","beta5.txt","beta6.txt",
                       "beta7.txt","trtAct.txt","trtPre.txt","alpha.txt",
                       "power.txt")
    spatialParNames = ("$\\theta_0$", "$\\theta_{1,0}$", "$\\theta_{1,1}$",
                       "$\\theta_{1,2}$", "$\\theta_{1,3}$", "$\\theta_{2,0}$",
                       "$\\theta_{2,1}$", "$\\theta_{2,2}$", "$\\theta_{2,3}$",
                       "$\\theta_3$", "$\\theta_4$", "$\\theta_5$",
                       "$\\theta_6$")


    ## get edge values
    edgePostDist = getPostDist(edgeSamplesFile, edgeSamplesOrder)

    edgeSimVals = getSimVals(edgeDir, edgeParFiles)

    assert len(edgePostDist) == len(edgeSimVals)
    assert len(edgePostDist) == len(edgeParNames)

    print "edge"
    for name, (postMean, postL95, postU95), simVal in zip(edgeParNames,
                                                          edgePostDist,
                                                          edgeSimVals):
        print "%s & %0.2f & [ %0.2f, %0.2f] & %0.2f" % (name, postMean, postL95,
                                                        postU95, simVal)


    ## get spatial values
    fns = {12: np.exp}
    spatialPostDist = getPostDist(spatialSamplesFile, spatialSamplesOrder, fns)

    spatialSimVals = getSimVals(spatialDir, spatialParFiles)
    for i in fns:
        spatialSimVals[i] = fns[i](spatialSimVals[i])

    assert len(spatialPostDist) == len(spatialSimVals)
    assert len(spatialPostDist) == len(spatialParNames)

    print "spatial"
    for name, (postMean, postL95, postU95), simVal in zip(spatialParNames,
                                                          spatialPostDist,
                                                          spatialSimVals):
        print "%s & %0.2f & [ %0.2f, %0.2f] & %0.2f" % (name, postMean, postL95,
                                                        postU95, simVal)


if __name__ == "__main__":
    main()
