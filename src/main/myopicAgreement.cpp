#include <gflags/gflags.h>
#include <glog/logging.h>
#include <iostream>
#include <iomanip>
#include "runM1MlesWns.hpp"

using namespace google;
using namespace gflags;


DEFINE_string(srcDir,"","Path to source directory");

int main(int argc, char ** argv){
    InitGoogleLogging(argv[0]);
    ParseCommandLineFlags(&argc,&argv,true);

    njm::sett.setup(std::string(argv[0]),FLAGS_srcDir);

    njm::toFile(false,
            njm::sett.datExt("edgeToEdge_flag_",".txt"));

    // typedef ModelTimeExpCavesGDistTrendPowCon MG;
    typedef Model2GravityEDist MG;

    typedef MG ME;

    typedef System<MG,ME> S;

    typedef MyopicAgent<ME> MA;

    typedef WnsFeatures3<ME> F;
    typedef RankAgent<F,ME> RA;

    // S s;
    S s("obsData.txt");
    s.setEdgeToEdge(false);
    s.modelGen_r.setType(MLES);
    s.modelEst_r = s.modelGen_r; // use truth

    int numReps = 100;

    MA ma;
    RA ra;

    ra.tp.jitterScale = -1;
    ra.setEdgeToEdge(false);

    ra.tp.putPar({1.,0.,0.});

    ma.setEdgeToEdge(false);

    std::vector<std::vector<int> > rankTrt;
    std::vector<std::vector<int> > myopicTrt;

    for (int i = 0; i < numReps; ++i) {
        njm::resetSeed(i);
        s.revert();
        for (int j = 0; j < s.fD.finalT; ++j) {
            if (j >= s.fD.trtStart) {
                // myopic
                std::fill(s.tD.a.begin(),s.tD.a.end(),0);
                std::fill(s.tD.p.begin(),s.tD.p.end(),0);
                ma.applyTrt(s.sD,s.tD,s.fD,s.dD,s.modelEst);
                s.updateStatus();

                myopicTrt.push_back(s.sD.status);

                // rank
                std::fill(s.tD.a.begin(),s.tD.a.end(),0);
                std::fill(s.tD.p.begin(),s.tD.p.end(),0);
                ra.applyTrt(s.sD,s.tD,s.fD,s.dD,s.modelEst);
                s.updateStatus();

                rankTrt.push_back(s.sD.status);
            }
            s.nextPoint();
        }
    }

    njm::toFile(njm::toString(rankTrt, "\n", ""),
            njm::sett.datExt("rank_trt_",".txt"));

    njm::toFile(njm::toString(myopicTrt, "\n", ""),
            njm::sett.datExt("myopic_trt_",".txt"));


    return 0;
}
