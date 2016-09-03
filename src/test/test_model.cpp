#include <gtest/gtest.h>
#include <glog/logging.h>
#include <boost/filesystem.hpp>
#include "settings.hpp"

TEST(TestModel,TestInit) {
}

TEST(TestModel,TestRead) {
}

TEST(TestModel,TestSave) {
}

TEST(TestModel,TestLinScale) {
}

TEST(TestModel,TestSetType) {
}

TEST(TestModel,TestSetFixSample) {
}

TEST(TestModel,TestGetType) {
}

TEST(TestModel,TestSetEdgeToEdge) {
}

TEST(TestModel,TestGetEdgeToEdge) {
}

TEST(TestModel,TestInfProbs) {
}

TEST(TestModel,TestRevProbs) {
}

TEST(TestModel,TestSetFill) {
}

TEST(TestModel,TestModFill) {
}

TEST(TestModel,TestSetQuick) {
}

TEST(TestModel,TestOneOnOne) {
}

TEST(TestModel,TestGetPar) {
}

TEST(TestModel,TestPutPar) {
}

TEST(TestModel,TestSetPar) {
}

TEST(TestModel,TestPartial) {
}

TEST(TestModel,TestPartial2) {
}

TEST(TestModel,TestSetFisher) {
}

TEST(TestModel,TestSample) {
}

TEST(TestModel,TestRevert) {
}

TEST(TestModel,TestFit) {
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  char fileName[32] = "test_model";
  boost::filesystem::path temp = boost::filesystem::unique_path();
  char srcDir[32];
  std::strcpy(srcDir,temp.native().c_str());

  char check[32] = "Y";

  char* pseudoArgv[3];
  pseudoArgv[0] = fileName;
  pseudoArgv[1] = srcDir;
  pseudoArgv[2] = check;
  njm::sett.set(3,pseudoArgv);
  return RUN_ALL_TESTS();
}
