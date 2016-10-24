#include <gflags/gflags.h>
#include <glog/logging.h>

#include "runM1Mles.hpp"

#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>

using namespace google;
using namespace gflags;

template <typename T>
class GradientChecker {
public:
    System<T,T> * system;
    T * m;
    int gradientVar;
};

template <typename T>
class HessianChecker {
public:
    System<T,T> * system;
    T * m;
    int gradientVar;
    int hessianVar;
};

template <typename T>
double f (double x, void * params) {
    GradientChecker<T> * gc = static_cast<GradientChecker<T>*>(params);
    std::vector<double> par = gc->m->getPar();
    par.at(gc->gradientVar) = x;
    gc->m->putPar(par.begin());
    return gc->m->logll(gc->system->sD,gc->system->tD,
            gc->system->fD,gc->system->dD);
}

template <typename T>
double fGrad (double x, void * params) {
    HessianChecker<T> * hc = static_cast<HessianChecker<T>*>(params);
    std::vector<double> par = hc->m->getPar();
    par.at(hc->hessianVar) = x;
    hc->m->putPar(par.begin());
    return hc->m->logllGrad(hc->system->sD,hc->system->tD,
            hc->system->fD,hc->system->dD).at(hc->gradientVar);
}


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
            typedef Model2EdgeToEdge MG;

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

            r = 6;
            njm::resetSeed(r);
            s.reset(starts[r]);

            for(t=s.sD.time; t<12; ++t) {
                if(t>=s.fD.trtStart && s.sD.numNotInfec > 0){
                    s.modelEst.fit(s.sD,s.tD,s.fD,s.dD,
                            t > s.fD.trtStart);

                    ma.applyTrt(s.sD,s.tD,s.fD,s.dD,
                            s.modelEst);
                }

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

            gsl_multimin_fdfminimizer_set(sfdf,&my_func,x,0.01,0.1);

            int iter = 0;
            int status;
            const int maxIter = 100;
            do{
                iter++;
                status = gsl_multimin_fdfminimizer_iterate(sfdf);

                if(status)
                    break;

                std::vector<double> currentX;
                for (pi = 0; pi < int(s.modelEst.numPars); ++pi) {
                    currentX.push_back(gsl_vector_get(sfdf->x,pi));
                }
                std::cout << std::endl
                          << "par: " << njm::toString(currentX," ","");

                status = gsl_multimin_test_gradient(sfdf->gradient,1e-6);

            }while(status == GSL_CONTINUE && iter < maxIter);

            std::vector<double> gradVals;
            for (pi = 0; pi < int(s.modelEst.numPars); ++pi) {
                gradVals.push_back(gsl_vector_get(sfdf->gradient,pi));
            }


            if (status != GSL_CONTINUE || status != GSL_SUCCESS) {
                LOG(ERROR)
                    << std::endl
                    << "status: " << status << std::endl
                    << "iter: " << iter << std::endl
                    << "seed: " << njm::getSeed() << std::endl
                    << "numInfected: " << s.sD.numInfected << std::endl
                    << "time: " << s.sD.time << std::endl
                    << "gradient check: "
                    << gsl_multimin_test_gradient(sfdf->gradient,0.1)
                    << std::endl
                    << "f: " << sfdf->f << std::endl
                    << "gradient: " << njm::toString(gradVals," ","")
                    << std::endl
                    << "starting: " << njm::toString(startingVals," ","")
                    << std::endl;
            }


            /// check gradient
            std::vector<double> par = s.modelEst.getPar();

            const double val = s.modelEst.logll(s.sD,s.tD,s.fD,s.dD);

            const std::vector<double> gradVal = s.modelEst.logllGrad(
                    s.sD,s.tD,s.fD,s.dD);

            const std::pair<double,std::vector<double> > both =
                s.modelEst.logllBoth(s.sD,s.tD,s.fD,s.dD);

            const double eps = 1e-6;

            CHECK_EQ(val,both.first);

            for (int i = 0; i < s.modelEst.numPars; ++i) {
                LOG(INFO) << "Checking par " << i;
                s.modelEst.putPar(par.begin());

                GradientChecker<ME> gc;
                gc.system = &s;
                gc.m = &s.modelEst;
                gc.gradientVar = i;

                gsl_function F;
                F.function = &f<ME>;
                F.params = &gc;

                double result;
                double abserr;
                gsl_deriv_central(&F,par.at(i),1e-8,&result,&abserr);

                CHECK_EQ(gradVal.at(i),both.second.at(i));
                CHECK_LT(std::abs(result-gradVal.at(i)),eps)
                    << "index: " << i << std::endl
                    << "par: " << njm::toString(par," ","") << std::endl
                    << "result: " << result << std::endl
                    << "gradVal: " << gradVal.at(i) << std::endl;
            }


        } else {
            LOG(FATAL) << "Supposed to debug edge-to-edge spread.";

        }

        njm::sett.clean();
    }
    return 0;
}
