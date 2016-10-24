#include <gflags/gflags.h>
#include <glog/logging.h>
#include "runM1Mles.hpp"

using namespace google;
using namespace gflags;


DEFINE_string(srcDir,"","Path to source directory");
DEFINE_bool(edgeToEdge,false,"Edge to edge transmission");
DEFINE_bool(dryRun,false,"Do not execute main");

int main(int argc, char ** argv){
    InitGoogleLogging(argv[0]);
    ParseCommandLineFlags(&argc,&argv,true);
    if(!FLAGS_dryRun) {
        njm::sett.setup(std::string(argv[0]),FLAGS_srcDir);

        njm::toFile(FLAGS_edgeToEdge,
                njm::sett.datExt("edgeToEdge_flag_",".txt"));

        if(FLAGS_edgeToEdge) {
            LOG(FATAL) << "Supposed to debug spatial spread.";
        } else {
            // typedef ModelTimeExpCavesGPowGDistTrendPowCon MG;
            typedef Model2GravityEDist MG;

            typedef MG ME;

            typedef System<MG,ME> S;

            typedef MyopicAgent<ME> MA;

            typedef FitOnlyRunner<S,MA> R_MA;


            S s;
            s.setEdgeToEdge(FLAGS_edgeToEdge);
            s.modelGen_r.setType(MLES);
            s.modelEst_r.setType(MLES);

            int numReps = 100;
            Starts starts(numReps,s.fD.numNodes);

            MA ma;

            ma.setEdgeToEdge(FLAGS_edgeToEdge);

            R_MA r_ma;

            int r,t;

            r = 65;
            njm::resetSeed(r);
            s.reset(starts[r]);

            for(t=s.sD.time; t<s.fD.trtStart; ++t) {
                s.updateStatus();

                s.nextPoint();
            }

            s.modelEst.fit(s.sD,s.tD,s.fD,s.dD,
                    t > s.fD.trtStart);

        }

        njm::sett.clean();
    }
    return 0;
}
