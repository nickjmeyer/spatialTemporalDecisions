#include <gtest/gtest.h>
#include <glog/logging.h>
#include <boost/filesystem.hpp>
#include "settings.hpp"
#include "data.hpp"
#include "system.hpp"
#include "paramIntercept.hpp"
#include "model.hpp"

class TestModel : public ::testing::Test {
public:
  ModelBase * m;

  SimData * sD;
  TrtData * tD;
  FixedData * fD;
  DynamicData * dD;

  TestModel() {


    // m = new ModelBase("TestModel",{new ParamIntercept},*fD);
  }

  ~TestModel() {
    // delete m;
    // delete sD;
    // delete tD;
    // delete fD;
    // delete dD;
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


  // setup the fake network
  std::vector<int> fips;
  fips.push_back(0);
  fips.push_back(1);
  fips.push_back(2);
  fips.push_back(3);
  njm::toFile(fips,njm::sett.srcExt("fips.txt"));

  std::vector<float> gDist;
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



  return RUN_ALL_TESTS();
}
