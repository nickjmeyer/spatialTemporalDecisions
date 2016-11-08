#include <gflags/gflags.h>
#include <gtest/gtest.h>
#include <glog/logging.h>
#include <boost/filesystem.hpp>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <algorithm>
#include "runStats.hpp"
#include "utilities.hpp"

#include <random>

using namespace google;
using namespace gflags;


const float eps = 1e-12;


TEST(TestRunStats,TestSampleValues) {
    std::vector<double> v = {1,1.1,2.3,4.2,3.6,7.8,2.1,7.2,-2,3};
    RunStats rs;

    double sum = 0.0;
    for (int i = 0;  i < v.size(); ++i) {
        rs.update(v.at(i));
        sum += v.at(i);
    }
    EXPECT_NEAR(rs.sMean(),sum/v.size(),eps) << "Mean not correct";
    EXPECT_NEAR(rs.sVar(),njm::sampVar(v),eps) << "Var not correct";
    EXPECT_NEAR(rs.sSd(),std::sqrt(njm::sampVar(v)),eps) << "Sd not correct";
    EXPECT_NEAR(rs.seMean(),std::sqrt(njm::sampVar(v)/v.size()),eps)
        << "Se not correct";
}


TEST(TestRunStats,TestReservoir) {
    std::vector<double> value;
    std::vector<unsigned int> index;

    const unsigned int n = 10000;

    std::mt19937 rng(0);
    std::uniform_real_distribution<double> dis(0,1);
    for (unsigned int i = 0; i < n; ++i) {
        value.push_back(dis(rng));
        index.push_back(i);
    }

    std::random_shuffle(index.begin(),index.end());

    RunStats rs_no_res_iter;
    RunStats rs_no_res_all;
    RunStats rs_res_iter;
    RunStats rs_res_all;

    for (int i = 0; i < n; i++) {
        rs_no_res_iter.update(value.at(i));
        rs_res_iter.update(index.at(i),value.at(i));
    }
    rs_no_res_all.update(value);
    rs_res_all.update(index,value);

    rs_res_iter.squash();
    rs_res_all.squash();

    // mean
    CHECK_EQ(rs_no_res_all.sMean(),rs_no_res_iter.sMean());
    CHECK_EQ(rs_res_all.sMean(),rs_res_iter.sMean());
    CHECK_NEAR(rs_no_res_all.sMean(),rs_res_all.sMean(),eps);
    CHECK_NEAR(rs_no_res_all.sMean(),0.5,
            1./std::sqrt(static_cast<double>(n)));

    // var
    CHECK_EQ(rs_no_res_all.sVar(),rs_no_res_iter.sVar());
    CHECK_EQ(rs_res_all.sVar(),rs_res_iter.sVar());
    CHECK_NEAR(rs_no_res_all.sVar(),rs_res_all.sVar(),eps);
    CHECK_NEAR(rs_no_res_all.sVar(),1./12.,
            1./std::sqrt(static_cast<double>(n)));

    // sd
    CHECK_EQ(rs_no_res_all.sSd(),rs_no_res_iter.sSd());
    CHECK_EQ(rs_res_all.sSd(),rs_res_iter.sSd());
    CHECK_NEAR(rs_no_res_all.sSd(),rs_res_all.sSd(),eps);
    CHECK_NEAR(rs_no_res_all.sSd(),std::sqrt(1./12.),
            1./std::sqrt(static_cast<double>(n)));

    // se
    CHECK_EQ(rs_no_res_all.seMean(),rs_no_res_iter.seMean());
    CHECK_EQ(rs_res_all.seMean(),rs_res_iter.seMean());
    CHECK_NEAR(rs_no_res_all.seMean(),rs_res_all.seMean(),eps);
    CHECK_NEAR(rs_no_res_all.seMean(),std::sqrt(1./12./n),1./n);

    // run again to see if squash resets properly
    for (int i = 0; i < n; i++) {
        rs_no_res_iter.update(value.at(i));
        rs_res_iter.update(index.at(i),value.at(i));
    }
    rs_no_res_all.update(value);
    rs_res_all.update(index,value);

    rs_res_iter.squash();
    rs_res_all.squash();

    // mean
    CHECK_EQ(rs_no_res_all.sMean(),rs_no_res_iter.sMean());
    CHECK_EQ(rs_res_all.sMean(),rs_res_iter.sMean());
    CHECK_NEAR(rs_no_res_all.sMean(),rs_res_all.sMean(),eps);
    CHECK_NEAR(rs_no_res_all.sMean(),0.5,1./std::sqrt(2.*n));

    // var
    CHECK_EQ(rs_no_res_all.sVar(),rs_no_res_iter.sVar());
    CHECK_EQ(rs_res_all.sVar(),rs_res_iter.sVar());
    CHECK_NEAR(rs_no_res_all.sVar(),rs_res_all.sVar(),eps);
    CHECK_NEAR(rs_no_res_all.sVar(),1./12.,1./std::sqrt(2.*n));

    // sd
    CHECK_EQ(rs_no_res_all.sSd(),rs_no_res_iter.sSd());
    CHECK_EQ(rs_res_all.sSd(),rs_res_iter.sSd());
    CHECK_NEAR(rs_no_res_all.sSd(),rs_res_all.sSd(),eps);
    CHECK_NEAR(rs_no_res_all.sSd(),std::sqrt(1./12.),1./std::sqrt(2.*n));

    // se
    CHECK_EQ(rs_no_res_all.seMean(),rs_no_res_iter.seMean());
    CHECK_EQ(rs_res_all.seMean(),rs_res_iter.seMean());
    CHECK_NEAR(rs_no_res_all.seMean(),rs_res_all.seMean(),eps);
    CHECK_NEAR(rs_no_res_all.seMean(),std::sqrt(1./12./(2*n)),1./(2.*n));
}


int main(int argc, char **argv) {
    ParseCommandLineFlags(&argc,&argv,true);
    ::testing::InitGoogleTest(&argc, argv);

    int ret = RUN_ALL_TESTS();
    return ret;
}
