#include <gflags/gflags.h>
#include <glog/logging.h>
#include "tuneGen.hpp"

using namespace google;
using namespace gflags;


DEFINE_string(srcDir,"","Path to source directory");
DEFINE_bool(edgeToEdge,false,"Edge to edge transmission");
DEFINE_string(outDir,"","Path to save parameters");
DEFINE_bool(dryRun,false,"Do not execute main");

template <typename T>
void copyParams(const boost::filesystem::path path) {
    System<T,T> s;
    s.modelGen_r.read();
    s.modelGen_r.save_to(path);
}

// double getDPow(const double & power, const double & alpha,
//   const std::vector<double> & caves){
//   double meanCaves = std::accumulate(caves.begin(),caves.end(),0);
//   meanCaves /= double(caves.size());

//   double dPow = std::log(2.0)*std::pow(meanCaves,2.0*std::exp(power))/alpha;
//   dPow += 1.0;
//   dPow = std::log(dPow);
//   dPow /= std::log(2.0);

//   return(dPow);
// }





template <class S, class NT,class RN, class MG>
double TuneGenNT(S s, const int numReps, const Starts & starts){
    NT nt;
    RN rn;

    double goal = 0.7;
    njm::message("Goal: " + njm::toString(goal,""));

    int numYears = s.fD.finalT;
    double tol = 1e-3;

    // s.modelGen_r.setPar("gPow",getDPow(s.modelGen_r.getPar({"power"})[0],
    //     s.modelGen_r.getPar({"alpha"})[0],
    //     s.fD.caves));

    std::vector<double> par = s.modelGen_r.getPar();
    s.modelEst_r.putPar(par.begin());
    s.revert();


    double val = rn.run(s,nt,numReps,numYears,starts).sMean();
    double scale = 1.025, shrink = .9;
    int above = int(val > goal);
    int iter = 0;

    printf("Iter: %05d  >>>  Current value: %08.6f\r",
            ++iter, val);

    while(std::abs(val - goal) > tol){
        if(val > goal){
            if(!above)
                scale*=shrink;

            // s.modelGen_r.linScale(1.0 + scale);

            std::vector<double> curIntcp = s.modelGen_r.getPar({"intcp"});
            CHECK_EQ(curIntcp.size(),1) << "more than one intercept was returned";
            curIntcp.at(0) -= scale;
            s.modelGen_r.setPar("intcp",curIntcp.at(0));

            above = 1;
        }
        else{
            if(above)
                scale*=shrink;

            std::vector<double> curIntcp = s.modelGen_r.getPar({"intcp"});
            CHECK_EQ(curIntcp.size(),1) << "more than one intercept was returned";
            curIntcp.at(0) += scale;
            s.modelGen_r.setPar("intcp",curIntcp.at(0));
            // s.modelGen_r.linScale(1.0/(1.0 + scale));

            above = 0;
        }

        par = s.modelGen_r.getPar();
        s.modelEst_r.putPar(par.begin());
        s.revert();


        // s.modelGen_r.setPar("gPow",getDPow(s.modelGen_r.getPar({"power"})[0],
        //     s.modelGen_r.getPar({"alpha"})[0],
        //     s.fD.caves));



        val = rn.run(s,nt,numReps,numYears,starts).sMean();
        printf("Iter: %05d  >>>  Current value: %08.6f  (%08.6f)\r",
                ++iter, val, scale);
        // std::cout << std::endl
        //           << njm::toString(s.modelGen_r.getPar()," ","") << std::endl;
        fflush(stdout);
    }

    s.modelGen_r.save();

    njm::message("Est. goal: " + njm::toString(val,""));


    njm::message("par: " + njm::toString(s.modelGen_r.getPar()," ",""));

    return(val);
}


template <class S, class MA, class RM, class NT, class RN>
double TuneGenMA(S s, const int numReps, const Starts & starts){
    NT nt;
    RN rn;

    MA ma;
    RM rm;

    int numYears = s.fD.finalT;

    double atTrtStart = rn.run(s,nt,numReps,s.fD.trtStart,starts).sMean();
    double atFinalT = rn.run(s,nt,numReps,numYears,starts).sMean();

    double goal = atTrtStart + 0.05*(atFinalT - atTrtStart);
    njm::message("Goal: " + njm::toString(goal,""));
    double tol = 1e-3;

    std::vector<double> par;
    double trt = 1.0;

    s.modelGen_r.setPar(std::vector<std::string>({"trtAct","trtPre"}),trt);
    par = s.modelGen_r.getPar();
    s.modelEst_r.putPar(par.begin());
    s.revert();

    double val = rm.run(s,ma,numReps,numYears,starts).sMean();
    double scale = 1.1, shrink = .9;
    int above = int(val > goal);
    int iter = 0;


    printf("Iter: %05d  >>>  Curr value: %08.6f  ===  Curr Trt: %08.6f\r",
            ++iter, val, trt);

    while(std::abs(val - goal) > tol){
        if(val > goal){
            if(!above)
                scale*=shrink;

            trt *= 1.0 + scale;

            above = 1;
        }
        else{
            if(above)
                scale*=shrink;

            trt *= 1.0/(1.0 + scale);

            above = 0;
        }


        s.modelGen_r.setPar(std::vector<std::string>({"trtAct","trtPre"}),trt);
        par = s.modelGen_r.getPar();
        s.modelEst_r.putPar(par.begin());
        s.revert();

        // std::cout << "par: " << njm::toString(par," ","\n");
        // par = s.modelGen.getPar({"trtAct","trtPre"});
        // std::cout << "par: " << njm::toString(par," ","\n");


        val = rm.run(s,ma,numReps,numYears,starts).sMean();
        printf("Iter: %05d  >>>  Curr value: %08.6f  ===  Curr Trt: %08.6f\r",
                ++iter, val, trt);
        fflush(stdout);
    }

    s.modelGen_r.save();

    double priorMeanTrt = (s.modelGen_r.getPar({"trtAct"})[0]
            + s.modelGen_r.getPar({"trtPre"})[0])/2.0;
    priorMeanTrt *= 4.0;
    njm::toFile(priorMeanTrt,njm::sett.srcExt("priorTrtMean.txt"),
            std::ios_base::out);


    njm::message("Est. goal: " + njm::toString(val,""));

    njm::message("par: " + njm::toString(s.modelGen_r.getPar()," ",""));

    return(val);
}


// template <class S, class PA, class RP>
// double TuneGenPA(S & s,const int numReps, const Starts & starts){
//   double trtSize = s.modelGen.tuneTrt(s.fD);

//   putActTrt(trtSize,s.modelGen_r,s.fD);
//   putPreTrt(trtSize,s.modelGen_r,s.fD);
//   putActTrt(trtSize,s.modelEst_r,s.fD);
//   putPreTrt(trtSize,s.modelEst_r,s.fD);

//   PA pa;
//   RP rp;

//   return rp.run(s,pa,numReps,s.fD.finalT,starts).sMean();
// }


int main(int argc, char ** argv){
    InitGoogleLogging(argv[0]);
    ParseCommandLineFlags(&argc,&argv,true);
    if(!FLAGS_dryRun) {
        njm::sett.setup(std::string(argv[0]),FLAGS_srcDir);

        const boost::filesystem::path path(FLAGS_outDir);

        if(FLAGS_edgeToEdge) {
            copyParams<Model2EdgeToEdge>(path);

            // typedef ModelTimeExpCavesGPowGDistTrendPowCon MG;

            typedef Model2EdgeToEdge MG;
            typedef MG ME;

            typedef System<MG,ME> S;
            typedef NoTrt<ME> NT;
            typedef ProximalAgent<ME> PA;
            typedef MyopicAgent<ME> MA;

            typedef AllAgent<ME> AA;

            typedef ToyFeatures5<ME> F;
            typedef RankAgent<F,ME> RA;

            typedef VanillaRunnerNS<S,NT> RN;
            typedef VanillaRunnerNS<S,PA> RP;
            typedef VanillaRunnerNS<S,MA> RM;
            typedef VanillaRunnerNS<S,RA> RR;

            typedef VanillaRunnerNS<S,AA> R_AA;

            S s;
            s.setEdgeToEdge(FLAGS_edgeToEdge);
            s.modelEst_r = s.modelGen_r;
            s.revert();

            int numReps = 500;
            Starts starts(numReps,s.fD.numNodes);

            MA ma;
            PA pa;
            RP rp;

            RA ra;
            RM rm;
            RR rr;

            pa.setEdgeToEdge(FLAGS_edgeToEdge);
            ra.setEdgeToEdge(FLAGS_edgeToEdge);
            ma.setEdgeToEdge(FLAGS_edgeToEdge);
            // ra.reset();

            ra.tp.jitterScale = -1.;

            njm::message("Tuning Intercept");

            double valNT = TuneGenNT<S,NT,RN,MG>(s,numReps,starts);
            s.modelGen_r.read();
            s.modelEst_r.read();
            s.revert();

            njm::message("Tuning Treatment");

            double valAA = TuneGenMA<S,AA,R_AA,NT,RN>(s,numReps,starts);
            s.modelGen_r.read();
            s.modelEst_r.read();
            s.revert();

            double valMA = rm.run(s,ma,numReps,s.fD.finalT,starts).sMean();

            double valPA = rp.run(s,pa,numReps,s.fD.finalT,starts).sMean();

            double valRA = rr.run(s,ra,numReps,s.fD.finalT,starts).sMean();

            njm::message(" valNT: " + njm::toString(valNT,"") +
                    "\n" +
                    " valPA: " + njm::toString(valPA,"") +
                    "\n" +
                    " valMA: " + njm::toString(valMA,"") +
                    "\n" +
                    " valRA: " + njm::toString(valRA,"") +
                    "\n" +
                    " valAA: " + njm::toString(valAA,""));

        } else {
            copyParams<Model2GravityEDist>(path);

            // typedef ModelTimeExpCavesGPowGDistTrendPowCon MG;

            typedef Model2GravityEDist MG;
            typedef MG ME;

            typedef System<MG,ME> S;
            typedef NoTrt<ME> NT;
            typedef ProximalAgent<ME> PA;
            typedef MyopicAgent<ME> MA;

            typedef AllAgent<ME> AA;

            typedef ToyFeatures5<ME> F;
            typedef RankAgent<F,ME> RA;

            typedef VanillaRunnerNS<S,NT> RN;
            typedef VanillaRunnerNS<S,PA> RP;
            typedef VanillaRunnerNS<S,MA> RM;
            typedef VanillaRunnerNS<S,RA> RR;

            typedef VanillaRunnerNS<S,AA> R_AA;

            S s;
            s.setEdgeToEdge(FLAGS_edgeToEdge);
            s.modelEst_r = s.modelGen_r;
            s.revert();

            int numReps = 500;
            Starts starts(numReps,s.fD.numNodes);

            MA ma;
            PA pa;
            RP rp;

            RA ra;
            RM rm;
            RR rr;

            pa.setEdgeToEdge(FLAGS_edgeToEdge);
            ra.setEdgeToEdge(FLAGS_edgeToEdge);
            ma.setEdgeToEdge(FLAGS_edgeToEdge);
            // ra.reset();

            ra.tp.jitterScale = -1.;

            njm::message("Tuning Intercept");

            double valNT = TuneGenNT<S,NT,RN,MG>(s,numReps,starts);
            s.modelGen_r.read();
            s.modelEst_r.read();
            s.revert();

            njm::message("Tuning Treatment");

            double valAA = TuneGenMA<S,AA,R_AA,NT,RN>(s,numReps,starts);
            s.modelGen_r.read();
            s.modelEst_r.read();
            s.revert();

            double valMA = rm.run(s,ma,numReps,s.fD.finalT,starts).sMean();

            double valPA = rp.run(s,pa,numReps,s.fD.finalT,starts).sMean();

            double valRA = rr.run(s,ra,numReps,s.fD.finalT,starts).sMean();

            njm::message(" valNT: " + njm::toString(valNT,"") +
                    "\n" +
                    " valPA: " + njm::toString(valPA,"") +
                    "\n" +
                    " valMA: " + njm::toString(valMA,"") +
                    "\n" +
                    " valRA: " + njm::toString(valRA,"") +
                    "\n" +
                    " valAA: " + njm::toString(valAA,""));
        }

        njm::sett.clean();

    }

    return 0;
}
