#include <gflags/gflags.h>
#include <gtest/gtest.h>
#include <glog/logging.h>
#include <boost/filesystem.hpp>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include "settings.hpp"
#include "data.hpp"
#include "system.hpp"
#include "toyFeatures5.hpp"
#include "wnsFeatures3.hpp"
#include "paramIntercept.hpp"
#include "modelIntercept.hpp"

const float eps = 1e-6;

template <typename T>
class TestFeatures : public ::testing::Test {
public:
    T f0;
    T f1;

    ModelIntercept * m;

    System<ModelIntercept,ModelIntercept> * system;

    TestFeatures() {
        system = new System<ModelIntercept,ModelIntercept>();
        system->reset({0});
        std::vector<double> mockProbs(system->sD.numNotInfec,0.0);
        mockProbs.at(0) = 1.0;
        mockProbs.at(system->sD.numNotInfec-1) = 1.0;

        system->nextPoint(mockProbs);

        m = new ModelIntercept(system->fD);
        m->read();
    }

    ~TestFeatures() {
        delete m;
        delete system;
    }
};


typedef ::testing::Types<
    ToyFeatures5<ModelIntercept>,
    WnsFeatures3<ModelIntercept>
    > MyTypes;

TYPED_TEST_CASE(TestFeatures,MyTypes);

TYPED_TEST(TestFeatures,TestUpdateEdgeToEdge) {
    this->system->setEdgeToEdge(true);
    this->m->setEdgeToEdge(true);
    this->f0.tp.setEdgeToEdge(true);
    this->f1.tp.setEdgeToEdge(true);

    // precomp data
    this->f0.preCompData(this->system->sD, this->system->tD, this->system->fD,
            this->system->dD,*this->m);
    this->f1.preCompData(this->system->sD, this->system->tD, this->system->fD,
            this->system->dD,*this->m);

    // get f0
    this->f0.getFeatures(this->system->sD,this->system->tD,this->system->fD,
            this->system->dD,*this->m);

    // apply trt
    this->system->tD.a.at(this->system->sD.infected.at(0)) = 1;

    // update f0
    this->f0.updateFeatures(this->system->sD,this->system->tD,this->system->fD,
            this->system->dD,*this->m);
    // get f1
    this->f1.getFeatures(this->system->sD,this->system->tD,this->system->fD,
            this->system->dD,*this->m);

    {
        const arma::mat diff = arma::abs(this->f0.infFeat - this->f1.infFeat);
        arma::mat::const_iterator it;
        const arma::mat::const_iterator end = diff.end();
        int i;
        for (it = diff.begin(), i = 0; it != end; ++it,++i) {
            EXPECT_LT(*it,eps) << "element " << i << " is not equal";
        }
    }

    {
        const arma::mat diff = arma::abs(this->f0.notFeat - this->f1.notFeat);
        arma::mat::const_iterator it;
        const arma::mat::const_iterator end = diff.end();
        int i;
        for (it = diff.begin(), i = 0; it != end; ++it,++i) {
            EXPECT_LT(*it,eps) << "element " << i << " for not infected is not equal";
        }
    }
}


TYPED_TEST(TestFeatures,TestUpdateSpatial) {
    this->system->setEdgeToEdge(false);
    this->m->setEdgeToEdge(false);
    this->f0.tp.setEdgeToEdge(false);
    this->f1.tp.setEdgeToEdge(false);

    // precomp data
    this->f0.preCompData(this->system->sD, this->system->tD, this->system->fD,
            this->system->dD,*this->m);
    this->f1.preCompData(this->system->sD, this->system->tD, this->system->fD,
            this->system->dD,*this->m);

    // get f0
    this->f0.getFeatures(this->system->sD,this->system->tD,this->system->fD,
            this->system->dD,*this->m);

    // apply trt
    this->system->tD.a.at(this->system->sD.infected.at(0)) = 1;

    // update f0
    this->f0.updateFeatures(this->system->sD,this->system->tD,this->system->fD,
            this->system->dD,*this->m);
    // get f1
    this->f1.getFeatures(this->system->sD,this->system->tD,this->system->fD,
            this->system->dD,*this->m);

    {
        const arma::mat diff = arma::abs(this->f0.infFeat - this->f1.infFeat);
        arma::mat::const_iterator it;
        const arma::mat::const_iterator end = diff.end();
        int i;
        for (it = diff.begin(), i = 0; it != end; ++it,++i) {
            EXPECT_LT(*it,eps) << "element " << i << " is not equal";
        }
    }

    {
        const arma::mat diff = arma::abs(this->f0.notFeat - this->f1.notFeat);
        arma::mat::const_iterator it;
        const arma::mat::const_iterator end = diff.end();
        int i;
        for (it = diff.begin(), i = 0; it != end; ++it,++i) {
            EXPECT_LT(*it,eps) << "element " << i << " for not infected is not equal";
        }
    }
}


void fakeNetworkSetup() {
    // setup the fake network
    const int numNodes = 9;
    const int numRow = 3;
    const int numCol = 3;
    ASSERT_EQ(numRow*numCol,numNodes);

    std::vector<int> fips;
    for (int i = 0; i < numNodes; ++i) {
        fips.push_back(i);
    }
    njm::toFile(fips,njm::sett.srcExt("fips.txt"));

    std::vector<double> gDist,eDist;
    for (int rowI = 0; rowI < numRow; ++rowI) {
        for (int colI = 0; colI < numCol; ++colI) {
            for (int rowJ = 0; rowJ < numRow; ++rowJ) {
                for (int colJ = 0; colJ < numCol; ++colJ) {
                    const int rowDiff = rowI - rowJ;
                    const int colDiff = colI - colJ;
                    gDist.push_back(std::abs(rowDiff) + std::abs(colDiff));
                    eDist.push_back(
                            std::sqrt((double)(rowDiff * rowDiff + colDiff * colDiff)));
                }
            }
        }
    }
    njm::toFile(eDist,njm::sett.srcExt("eDist.txt"));
    njm::toFile(gDist,njm::sett.srcExt("gDist.txt"));

    std::vector<int> caves;
    for (int i = 0; i < numNodes; ++i) {
        caves.push_back(i+1);
    }
    njm::toFile(caves,njm::sett.srcExt("caves.txt"));

    std::vector<double> xcov;
    const int numCov = 2;
    for (int i = 0; i < numNodes; ++i) {
        for (int j = 0; j < numCov; ++j) {
            xcov.push_back(i + 0.1 * j);
        }
    }
    njm::toFile(xcov,njm::sett.srcExt("xcov.txt"));

    std::vector<int> network;
    for (int i = 0; i < numNodes; ++i) {
        for (int j = 0; j < numNodes; ++j) {
            if(std::abs(gDist.at(i*numNodes + j) - 1.0)
                    < eps || i == j) {
                network.push_back(1);
            } else {
                network.push_back(0);
            }
        }
    }
    njm::toFile(network,njm::sett.srcExt("network.txt"));

    std::vector<double> centroidsLong;
    std::vector<double> centroidsLat;
    for (int row = 0; row < numRow; ++row) {
        for (int col = 0; col < numCol; ++col) {
            centroidsLat.push_back(row);
            centroidsLong.push_back(col);
        }
    }

    njm::toFile(centroidsLong,njm::sett.srcExt("centroidsLong.txt"));
    njm::toFile(centroidsLat,njm::sett.srcExt("centroidsLat.txt"));
    njm::toFile(centroidsLong,njm::sett.srcExt("centroidsMdsLong.txt"));
    njm::toFile(centroidsLat,njm::sett.srcExt("centroidsMdsLat.txt"));

    std::vector<double> subGraph;
    for (int i = 0; i < numNodes; ++i) {
        subGraph.push_back(i);
    }
    njm::toFile(subGraph,njm::sett.srcExt("subGraph.txt"));

    std::vector<double> betweenness;
    for (int i = 0; i < numNodes; ++i) {
        betweenness.push_back(numNodes - i - 1);
    }
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
    const std::vector<std::string> parNames = {"intcp","trtAct","trtPre"};

    boost::filesystem::create_directories(njm::sett.srcExt("Param"+modName));
    for (int j = 0; j < parNames.size(); ++j) {
        const std::string parName = parNames.at(j);
        njm::toFile(njm::runif(-0.5,0.5),njm::sett.srcExt("Param"+modName+
                        "/"+parName+".txt"));
    }
}





int main(int argc, char **argv) {
    ::google::ParseCommandLineFlags(&argc,&argv,true);
    ::testing::InitGoogleTest(&argc, argv);
    const std::string fileName = "test_features";
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
