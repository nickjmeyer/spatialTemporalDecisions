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

            ////// debug model fitting
            std::vector<double> startingVals(s.modelEst.numPars,0.2);
            startingVals.at(0) = -3.0;

            ModelBaseFitObj fitObj(&s.modelEst,s.sD,s.tD,s.fD,s.dD);

            const gsl_multimin_fdfminimizer_type * T;
            gsl_multimin_fdfminimizer *sfdf;

            gsl_vector * x;
            x = gsl_vector_alloc(s.modelEst.numPars);
            int pi;
            for(pi = 0; pi < int(s.modelEst.numPars); ++pi){
                gsl_vector_set(x,pi,startingVals.at(pi));
            }

            gsl_multimin_function_fdf my_func;
            my_func.n = s.modelEst.numPars;
            my_func.f = objFn;
            my_func.df = objFnGrad;
            my_func.fdf = objFnBoth;
            my_func.params = &fitObj;

            T = gsl_multimin_fdfminimizer_vector_bfgs2;
            sfdf = gsl_multimin_fdfminimizer_alloc(T,s.modelEst.numPars);

            gsl_multimin_fdfminimizer_set(sfdf,&my_func,x,0.1,0.05);

            int iter = 0;
            int status;
            const int maxIter = 100;
            do{
                iter++;
                status = gsl_multimin_fdfminimizer_iterate(sfdf);

                if(status)
                    break;

                status = gsl_multimin_test_gradient(sfdf->gradient,1e-6);

            }while(status == GSL_CONTINUE && iter < maxIter);

        }

        njm::sett.clean();
    }
    return 0;
}
