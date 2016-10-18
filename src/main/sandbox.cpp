#include <gflags/gflags.h>
#include <glog/logging.h>
#include "runM1Mles.hpp"

DEFINE_string(srcDir,"","Path to source directory");
DEFINE_bool(edgeToEdge,false,"Edge to edge transmission");
DEFINE_bool(dryRun,false,"Do not execute main");

int main(int argc, char ** argv){
    ::google::InitGoogleLogging(argv[0]);
    ::google::ParseCommandLineFlags(&argc,&argv,true);
    if(!FLAGS_dryRun) {
        njm::sett.setup(std::string(argv[0]),FLAGS_srcDir);

        if(FLAGS_edgeToEdge) {
            // typedef ModelTimeExpCavesGDistTrendPowCon MG;
            typedef Model2EdgeToEdge MG;

            typedef MG ME;

            typedef System<MG,ME> S;

            typedef WnsFeatures3<ME> F;
            typedef RankAgent<F,ME> RA;

            typedef VanillaRunnerNS<S,RA> RR;

            S s("obsData.txt");
            s.setEdgeToEdge(FLAGS_edgeToEdge);
            s.modelEst_r = s.modelGen_r;
            s.revert();

            int numReps = 8;
            Starts starts("startingLocations.txt");

            RA ra;

            ra.setEdgeToEdge(FLAGS_edgeToEdge);
            // ra.reset();

            int r = numReps;
            int t;

            njm::resetSeed(r);
            s.reset(starts[r]);
            for(t=s.sD.time; t<s.fD.finalT; t++){
                if(t>=s.fD.trtStart && s.sD.numNotInfec > 0)
                    ra.applyTrt(s.sD,s.tD,s.fD,s.dD,
                            s.modelEst);

                s.updateStatus();

                s.nextPoint();

            }
            s.value();

            njm::message("value: " + njm::toString(s.value(),"",32,28));

        } else {
            // typedef ModelTimeExpCavesGDistTrendPowCon MG;
            // typedef Model2GravityEDist MG;

            // typedef MG ME;

            // typedef System<MG,ME> S;
            // typedef NoTrt<ME> NT;
            // typedef ProximalAgent<ME> PA;
            // typedef MyopicAgent<ME> MA;

            // typedef WnsFeatures3<ME> F;
            // typedef RankAgent<F,ME> RA;

            // typedef VanillaRunnerNS<S,NT> RN;
            // typedef VanillaRunnerNS<S,PA> RP;
            // typedef VanillaRunnerNS<S,MA> RM;
            // typedef VanillaRunnerNS<S,RA> RR;

            // S s("obsData.txt");
            // s.setEdgeToEdge(FLAGS_edgeToEdge);
            // s.modelEst_r = s.modelGen_r;
            // s.revert();

            // int numReps = 7;
            // Starts starts("startingLocations.txt");

            // NT nt;
            // MA ma;
            // PA pa;
            // RP rp;

            // RN rn;
            // RA ra;
            // RM rm;
            // RR rr;

            // pa.setEdgeToEdge(FLAGS_edgeToEdge);
            // ra.setEdgeToEdge(FLAGS_edgeToEdge);
            // ma.setEdgeToEdge(FLAGS_edgeToEdge);
            // // ra.reset();

            // double valRA = rr.run(s,ra,numReps,s.fD.finalT,starts).sMean();

            // njm::message(" valNT: " + njm::toString(valNT,"",32,28) +
            //         "\n" +
            //         " valPA: " + njm::toString(valPA,"",32,28) +
            //         "\n" +
            //         " valMA: " + njm::toString(valMA,"",32,28) +
            //         "\n" +
            //         " valRA: " + njm::toString(valRA,"",32,28));
        }

        njm::sett.clean();
    }

    return 0;
}
