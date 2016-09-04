#include <gtest/gtest.h>
#include <glog/logging.h>
#include <boost/filesystem.hpp>
#include "settings.hpp"
#include "data.hpp"
#include "system.hpp"
#include "paramIntercept.hpp"
#include "modelGravityGDist.hpp"

class TestModel : public ::testing::Test {
public:
  ModelGravityGDist * m;

  System<ModelGravityGDist,ModelGravityGDist> * system;

  TestModel() {
    system = new System<ModelGravityGDist,ModelGravityGDist>();

    m = new ModelGravityGDist(system->fD);
  }

  ~TestModel() {
    delete m;
    delete system;
  }
};

TEST_F(TestModel,TestInit) {
}

TEST_F(TestModel,TestRead) {
}

TEST_F(TestModel,TestSave) {
}

TEST_F(TestModel,TestLinScale) {
}

TEST_F(TestModel,TestSetType) {
}

TEST_F(TestModel,TestSetFixSample) {
}

TEST_F(TestModel,TestGetType) {
}

TEST_F(TestModel,TestSetEdgeToEdge) {
}

TEST_F(TestModel,TestGetEdgeToEdge) {
}

TEST_F(TestModel,TestInfProbs) {
}

TEST_F(TestModel,TestRevProbs) {
}

TEST_F(TestModel,TestSetFill) {
}

TEST_F(TestModel,TestModFill) {
}

TEST_F(TestModel,TestSetQuick) {
}

TEST_F(TestModel,TestOneOnOne) {
}

TEST_F(TestModel,TestGetPar) {
}

TEST_F(TestModel,TestPutPar) {
}

TEST_F(TestModel,TestSetPar) {
}

TEST_F(TestModel,TestPartial) {
}

TEST_F(TestModel,TestPartial2) {
}

TEST_F(TestModel,TestSetFisher) {
}

TEST_F(TestModel,TestSample) {
}

TEST_F(TestModel,TestRevert) {
}

TEST_F(TestModel,TestFit) {
}


void fakeNetworkSetup() {
  // setup the fake network
  std::vector<int> fips;
  fips.push_back(0);
  fips.push_back(1);
  fips.push_back(2);
  fips.push_back(3);
  njm::toFile(fips,njm::sett.srcExt("fips.txt"));

  std::vector<double> gDist;
  gDist.push_back(0.);
  gDist.push_back(1.);
  gDist.push_back(1.);
  gDist.push_back(2.);

  gDist.push_back(1.);
  gDist.push_back(0.);
  gDist.push_back(2.);
  gDist.push_back(1.);

  gDist.push_back(1.);
  gDist.push_back(2.);
  gDist.push_back(0.);
  gDist.push_back(1.);

  gDist.push_back(2.);
  gDist.push_back(1.);
  gDist.push_back(1.);
  gDist.push_back(0.);
  njm::toFile(gDist,njm::sett.srcExt("gDist.txt"));

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

  boost::filesystem::create_directories(njm::sett.srcExt("ParamGravityGDist"));
  njm::toFile(1.0,njm::sett.srcExt("ParamGravityGDist/intcp.txt"));
  njm::toFile(0.1,njm::sett.srcExt("ParamGravityGDist/beta0.txt"));
  njm::toFile(-0.1,njm::sett.srcExt("ParamGravityGDist/beta1.txt"));
  njm::toFile(0.1,njm::sett.srcExt("ParamGravityGDist/alpha.txt"));
  njm::toFile(0.1,njm::sett.srcExt("ParamGravityGDist/power.txt"));
  njm::toFile(0.5,njm::sett.srcExt("ParamGravityGDist/trtAct.txt"));
  njm::toFile(0.25,njm::sett.srcExt("ParamGravityGDist/trtPre.txt"));
}





int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  char fileName[32] = "test_model";
  boost::filesystem::path tempModel = boost::filesystem::temp_directory_path();
  tempModel += "/%%%%-%%%%-%%%%-%%%%";
  boost::filesystem::path temp = boost::filesystem::unique_path(tempModel);
  char srcDir[32];
  std::strcpy(srcDir,temp.native().c_str());

  char check[32] = "Y";

  char* pseudoArgv[3];
  pseudoArgv[0] = fileName;
  pseudoArgv[1] = srcDir;
  pseudoArgv[2] = check;
  njm::sett.set(3,pseudoArgv);

  fakeNetworkSetup();

  int ret = RUN_ALL_TESTS();
  return ret;
}
