#include <gflags/gflags.h>
#include <gtest/gtest.h>
#include <glog/logging.h>
#include <boost/filesystem.hpp>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include "settings.hpp"
#include "data.hpp"
#include "system.hpp"
#include "proximalAgent.hpp"
#include "paramIntercept.hpp"

TEST(TestModel,TestModFill) {
    System<ModelIntercept,ModelIntercept> s;
    s.modelEst_r = s.modelGen_r;
    s.revert();

    ProximalAgent<ModelIntercept> pa;

    s.reset({0});
    for (int i = 0; i < 100; ++i) {
        // without treatment
        for (int j = 0; j < s.sD.numInfected; ++j) {
            const int node = s.sD.infected.at(j);
            if (s.tD.a.at(node) == 0) {
                CHECK_EQ(s.sD.status.at(node),2);
            } else {
                CHECK_EQ(s.sD.status.at(node),3);
            }
        }
        for (int j = 0; j < s.sD.numNotInfec; ++j) {
            const int node = s.sD.notInfec.at(j);
            if(s.tD.p.at(node) == 0) {
                CHECK_EQ(s.sD.status.at(node),0);
            } else {
                CHECK_EQ(s.sD.status.at(node),1);
            }
        }

        // with treatment
        pa.applyTrt(s.sD,s.tD,s.fD,s.dD,s.modelEst);
        s.updateStatus();
        for (int j = 0; j < s.sD.numInfected; ++j) {
            const int node = s.sD.infected.at(j);
            if (s.tD.a.at(node) == 0) {
                CHECK_EQ(s.sD.status.at(node),2);
            } else {
                CHECK_EQ(s.sD.status.at(node),3);
            }
        }
        for (int j = 0; j < s.sD.numNotInfec; ++j) {
            const int node = s.sD.notInfec.at(j);
            if(s.tD.p.at(node) == 0) {
                CHECK_EQ(s.sD.status.at(node),0);
            } else {
                CHECK_EQ(s.sD.status.at(node),1);
            }
        }

        s.nextPoint();
    }

}

void fakeNetworkSetup() {
    // setup the fake network
    std::vector<int> fips;
    fips.push_back(0);
    fips.push_back(1);
    fips.push_back(2);
    fips.push_back(3);
    njm::toFile(fips,njm::sett.srcExt("fips.txt"));

    std::vector<double> dist;
    dist.push_back(0.);
    dist.push_back(1.);
    dist.push_back(1.);
    dist.push_back(2.);

    dist.push_back(1.);
    dist.push_back(0.);
    dist.push_back(2.);
    dist.push_back(1.);

    dist.push_back(1.);
    dist.push_back(2.);
    dist.push_back(0.);
    dist.push_back(1.);

    dist.push_back(2.);
    dist.push_back(1.);
    dist.push_back(1.);
    dist.push_back(0.);
    njm::toFile(dist,njm::sett.srcExt("gDist.txt"));
    njm::toFile(dist,njm::sett.srcExt("eDist.txt"));

    std::vector<int> caves;
    caves.push_back(1);
    caves.push_back(2);
    caves.push_back(3);
    caves.push_back(4);
    njm::toFile(caves,njm::sett.srcExt("caves.txt"));

    std::vector<double> xcov;
    xcov.push_back(0.);
    xcov.push_back(0.1);

    xcov.push_back(1.);
    xcov.push_back(1.1);

    xcov.push_back(2.);
    xcov.push_back(2.1);

    xcov.push_back(3.);
    xcov.push_back(3.1);
    njm::toFile(xcov,njm::sett.srcExt("xcov.txt"));

    std::vector<int> network;
    network.push_back(1);
    network.push_back(1);
    network.push_back(1);
    network.push_back(0);

    network.push_back(1);
    network.push_back(1);
    network.push_back(0);
    network.push_back(1);

    network.push_back(1);
    network.push_back(0);
    network.push_back(1);
    network.push_back(1);

    network.push_back(0);
    network.push_back(1);
    network.push_back(1);
    network.push_back(1);
    njm::toFile(network,njm::sett.srcExt("network.txt"));

    std::vector<double> centroidsLong;
    centroidsLong.push_back(0.);
    centroidsLong.push_back(1.);
    centroidsLong.push_back(0.);
    centroidsLong.push_back(1.);

    std::vector<double> centroidsLat;
    centroidsLat.push_back(0.);
    centroidsLat.push_back(0.);
    centroidsLat.push_back(1.);
    centroidsLat.push_back(1.);

    njm::toFile(centroidsLong,njm::sett.srcExt("centroidsLong.txt"));
    njm::toFile(centroidsLat,njm::sett.srcExt("centroidsLat.txt"));
    njm::toFile(centroidsLong,njm::sett.srcExt("centroidsMdsLong.txt"));
    njm::toFile(centroidsLat,njm::sett.srcExt("centroidsMdsLat.txt"));

    std::vector<double> subGraph;
    subGraph.push_back(0.);
    subGraph.push_back(1.);
    subGraph.push_back(2.);
    subGraph.push_back(3.);

    njm::toFile(subGraph,njm::sett.srcExt("subGraph.txt"));

    std::vector<double> betweenness;
    betweenness.push_back(3.);
    betweenness.push_back(2.);
    betweenness.push_back(1.);
    betweenness.push_back(0.);

    njm::toFile(betweenness,njm::sett.srcExt("betweenness.txt"));

    double priorTrtMean = 1.0;
    njm::toFile(priorTrtMean,njm::sett.srcExt("priorTrtMean.txt"));

    int trtStart = 3;
    njm::toFile(trtStart,njm::sett.srcExt("trtStart.txt"));

    int period = 1;
    njm::toFile(period,njm::sett.srcExt("period.txt"));

    int finalT = 6;
    njm::toFile(finalT,njm::sett.srcExt("finalT.txt"));

    const std::string modName = "Intercept";
    boost::filesystem::create_directories(njm::sett.srcExt("Param"+modName));
    // some pars aren't used in every model, but its cheap just to
    // write all parameters for all models rather than listing them
    // separately for each model
    njm::toFile(njm::runif(-0.5,0.5),njm::sett.srcExt("Param"+modName+
                    "/intcp.txt"));
    njm::toFile(njm::runif(-0.5,0.5),njm::sett.srcExt("Param"+modName+
                    "/trtAct.txt"));
    njm::toFile(njm::runif(-0.5,0.5),njm::sett.srcExt("Param"+modName+
                    "/trtPre.txt"));
}





int main(int argc, char **argv) {
    ::google::ParseCommandLineFlags(&argc,&argv,true);
    ::testing::InitGoogleTest(&argc, argv);
    const std::string fileName = "test_model";
    boost::filesystem::path tempModel = boost::filesystem::temp_directory_path();
    tempModel += "/%%%%-%%%%-%%%%-%%%%";
    boost::filesystem::path temp = boost::filesystem::unique_path(tempModel);
    const std::string srcDir(temp.native());

    njm::sett.setup(fileName,srcDir);

    fakeNetworkSetup();

    int ret = RUN_ALL_TESTS();
    njm::sett.clean();
    boost::filesystem::remove_all(temp);
    return ret;
}
