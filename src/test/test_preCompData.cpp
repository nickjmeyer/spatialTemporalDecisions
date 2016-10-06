#include <gtest/gtest.h>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include "tuneGenWNS.hpp"

DEFINE_string(srcDir,"","Path to source directory");

TEST(TestPreCompData,TestExpDistWeight) {
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

  const int cutoff = int(((double)s.fD.numNodes)/
    std::log((double)s.fD.numNodes));
  double sumProp = 0.0;
  double sumAll = 0.0;
  for (int i = 0;  i < uniqWeight.size(); ++i) {
    if (i < cutoff) {
      sumProp += uniqWeight.at(i);
    }
    sumAll += uniqWeight.at(i);
  }

  EXPECT_NEAR(sumProp,0.8*sumAll,1e-8);
}

int main(int argc, char **argv) {
  gflags::ParseCommandLineFlags(&argc,&argv,true);
  ::testing::InitGoogleTest(&argc, argv);
  const std::string fileName = "test_preCompData";
  const std::string srcDir("./data/wns");

  njm::sett.setup(fileName,srcDir);

  int ret = RUN_ALL_TESTS();
  njm::sett.clean();
  return ret;
}
