#include <gflags/gflags.h>
#include <gtest/gtest.h>
#include <glog/logging.h>
#include <boost/filesystem.hpp>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include "utilities.hpp"

const float eps = 1e-6;


TEST(TestUtilities,TestSampVar) {
  std::vector<double> v = {1,2,3,4,5};
  const double answer = (2.*2. + 1.*1. + 0. + 1.*1. + 2.*2.)/4.;
  EXPECT_NEAR(njm::sampVar(v),answer,1e-10);
}


int main(int argc, char **argv) {
  gflags::ParseCommandLineFlags(&argc,&argv,true);
  ::testing::InitGoogleTest(&argc, argv);

  int ret = RUN_ALL_TESTS();
  return ret;
}
