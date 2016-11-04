#include <gtest/gtest.h>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <algorithm>
#include "tuneGenWNS.hpp"

using namespace google;
using namespace gflags;


#ifndef STRINGIFY
#define STRINGIFY(x) #x
#endif

#ifndef TOSTRING
#define TOSTRING(x) STRINGIFY(x)
#endif




DEFINE_string(srcDir,"","Path to source directory");

TEST(TestPrecompData,TestExpDistWeightGradient) {
    std::vector<double> dVals;
    for (int i = 0; i < 100; ++i) {
        dVals.push_back(njm::runif01());
    }
    std::sort(dVals.begin(),dVals.end());

    ExpDistData edd;
    edd.dist = dVals;
    edd.proportion = 0.8;
    edd.cutoff = 20;

    gsl_function F;
    F.function = &expDistEval;
    F.params = &edd;

    double result;
    double abserr;

    const std::vector<double> rootVals = {-5.0,-1.0,0.0,1.0,5.0,10.0,100.0};
    for (int i = 0; i < rootVals.size(); ++i) {
        gsl_deriv_central(&F,rootVals.at(i),1e-8,&result,&abserr);
        EXPECT_NEAR(expDistGrad(rootVals.at(i),&edd),result,1e-6)
            << "Failed for root " << rootVals.at(i);
    }

}

TEST(TestPreCompData,TestExpDistWeightValue) {
    const bool edgeToEdge = false;

    typedef Model2GravityEDist MG;

    typedef MG ME;

    typedef System<MG,ME> S;

    typedef WnsFeatures3<ME> F;
    typedef RankAgent<F,ME> RA;

    typedef VanillaRunnerNS<S,RA> RR;

    S s("obsData.txt");

    std::vector<double> uniqWeight;
    for (int i = 0; i < s.fD.numNodes; ++i) {
        for (int j = (i+1); j < s.fD.numNodes; ++j) {
            uniqWeight.push_back(s.fD.expDistWeight.at(i*s.fD.numNodes + j));
        }
    }
    std::sort(uniqWeight.begin(),uniqWeight.end());
    std::reverse(uniqWeight.begin(),uniqWeight.end());

    const int cutoff = int(((double)uniqWeight.size())/
            std::log((double)uniqWeight.size()));
    double sumProp = 0.0;
    double sumAll = 0.0;
    for (int i = 0;  i < uniqWeight.size(); ++i) {
        if (i < cutoff) {
            sumProp += uniqWeight.at(i);
        }
        sumAll += uniqWeight.at(i);
    }

    EXPECT_NEAR(sumProp/sumAll,0.8,1e-8);
}

TEST(TestPreCompData,TestExpDistWeightNear) {
    const bool edgeToEdge = false;

    typedef Model2GravityEDist MG;

    typedef MG ME;

    typedef System<MG,ME> S;

    typedef WnsFeatures3<ME> F;
    typedef RankAgent<F,ME> RA;

    typedef VanillaRunnerNS<S,RA> RR;

    S s;

    for (int i = 0; i < s.fD.numNodes; ++i) {
        const std::vector<int> & near = s.fD.expDistWeightNear.at(i);
        const int numNear = near.size();
        double maxNearDist = -1.0;
        for (int j = 0; j < numNear; ++j) {
            const double distVal = s.fD.eDist.at(i*s.fD.numNodes + near.at(j));
            if(distVal > maxNearDist) {
                maxNearDist = distVal;
            }
        }
        ASSERT_GT(maxNearDist,0.0);

        for (int j = 0; j < s.fD.numNodes; ++j) {
            if(std::find(near.begin(), near.end(), j) == near.end()
                    && i != j) {
                EXPECT_GE(s.fD.eDist.at(i*s.fD.numNodes + j), maxNearDist);
            }
        }
    }
}


int main(int argc, char **argv) {
    ParseCommandLineFlags(&argc,&argv,true);
    ::testing::InitGoogleTest(&argc, argv);
    const std::string fileName = "test_preCompData";
    const std::string srcDir(TOSTRING(PROJECT_ROOT_DIR) "/data/wns");

    njm::sett.setup(fileName,srcDir);

    int ret = RUN_ALL_TESTS();
    njm::sett.clean();
    return ret;
}
