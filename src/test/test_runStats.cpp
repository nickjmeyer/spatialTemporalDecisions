#include <gflags/gflags.h>
#include <gtest/gtest.h>
#include <glog/logging.h>
#include <boost/filesystem.hpp>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include "runStats.hpp"
#include "utilities.hpp"

const float eps = 1e-6;


TEST(TestRunStats,TestSampleValues) {
  std::vector<double> v = {1,1.1,2.3,4.2,3.6,7.8,2.1,7.2,-2,3};
  RunStats rs;

  double sum = 0.0;
  for (int i = 0;  i < v.size(); ++i) {
    rs(v.at(i));
    sum += v.at(i);
  }
  EXPECT_NEAR(rs.sMean(),sum/v.size(),1e-10) << "Mean not correct";
  EXPECT_NEAR(rs.sVar(),njm::sampVar(v),1e-10) << "Var not correct";
  EXPECT_NEAR(rs.sSd(),std::sqrt(njm::sampVar(v)),1e-10) << "Sd not correct";
  EXPECT_NEAR(rs.seMean(),std::sqrt(njm::sampVar(v)/v.size()),1e-10)
    << "Se not correct";
}


int main(int argc, char **argv) {
    :;google::ParseCommandLineFlags(&argc,&argv,true);
  ::testing::InitGoogleTest(&argc, argv);

  int ret = RUN_ALL_TESTS();
  return ret;
}
