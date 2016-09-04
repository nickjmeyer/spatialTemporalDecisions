#include <gtest/gtest.h>
#include <glog/logging.h>
#include <boost/filesystem.hpp>
#include "settings.hpp"
#include "data.hpp"
#include "system.hpp"
#include "paramIntercept.hpp"
#include "modelGravityGDist.hpp"

const float eps = 1e-8;

class TestModel : public ::testing::Test {
public:
  static const double intcp;
  static const double alpha;
  static const double trtAct;
  static const double trtPre;

  ModelGDist * m;

  System<ModelGDist,ModelGDist> * system;

  TestModel() {
    system = new System<ModelGDist,ModelGDist>();
    system->reset({0});
    std::vector<double> mockProbs = {1.0,0.0,0.0};
    system->nextPoint(mockProbs);

    m = new ModelGDist(system->fD);
  }

  ~TestModel() {
    delete m;
    delete system;
  }
};


const double TestModel::intcp = 1.0;
const double TestModel::alpha = 0.15;
const double TestModel::trtAct = 0.5;
const double TestModel::trtPre = 0.25;


TEST_F(TestModel,TestInit) {
  EXPECT_EQ(this->m->numPars,4);
  EXPECT_EQ(this->m->name,"GDist");
  for (int i = 0; i < this->m->numPars; ++i) {
    EXPECT_NEAR(this->m->getPar().at(i),0.0,eps);
  }
}

TEST_F(TestModel,TestRead) {
  this->m->read();
  std::vector<double> par = this->m->getPar();
  EXPECT_NEAR(par.at(0),this->intcp,eps);
  EXPECT_NEAR(par.at(1),this->alpha,eps);
  EXPECT_NEAR(par.at(2),this->trtAct,eps);
  EXPECT_NEAR(par.at(3),this->trtPre,eps);
}

// TEST_F(TestModel,TestSave) {
// }

TEST_F(TestModel,TestLinScale) {
  {
    // unchanged
    const double scale = 1.0;
    this->m->read();
    this->m->linScale(scale);
    std::vector<double> par = this->m->getPar();
    EXPECT_NEAR(par.at(0),scale*this->intcp,eps);
    EXPECT_NEAR(par.at(1),scale*this->alpha,eps);
    EXPECT_NEAR(par.at(2),scale*this->trtAct,eps);
    EXPECT_NEAR(par.at(3),scale*this->trtPre,eps);
  }

  {
    // double
    const double scale = 2.0;
    this->m->read();
    this->m->linScale(scale);
    std::vector<double> par = this->m->getPar();
    EXPECT_NEAR(par.at(0),scale*this->intcp,eps);
    EXPECT_NEAR(par.at(1),scale*this->alpha,eps);
    EXPECT_NEAR(par.at(2),scale*this->trtAct,eps);
    EXPECT_NEAR(par.at(3),scale*this->trtPre,eps);
  }

  {
    // half
    const double scale = 0.5;
    this->m->read();
    this->m->linScale(scale);
    std::vector<double> par = this->m->getPar();
    EXPECT_NEAR(par.at(0),scale*this->intcp,eps);
    EXPECT_NEAR(par.at(1),scale*this->alpha,eps);
    EXPECT_NEAR(par.at(2),scale*this->trtAct,eps);
    EXPECT_NEAR(par.at(3),scale*this->trtPre,eps);
  }

  {
    // negative
    const double scale = -1.0;
    this->m->read();
    this->m->linScale(scale);
    std::vector<double> par = this->m->getPar();
    EXPECT_NEAR(par.at(0),scale*this->intcp,eps);
    EXPECT_NEAR(par.at(1),scale*this->alpha,eps);
    EXPECT_NEAR(par.at(2),scale*this->trtAct,eps);
    EXPECT_NEAR(par.at(3),scale*this->trtPre,eps);
  }

  {
    // negative double
    const double scale = -2.0;
    this->m->read();
    this->m->linScale(scale);
    std::vector<double> par = this->m->getPar();
    EXPECT_NEAR(par.at(0),scale*this->intcp,eps);
    EXPECT_NEAR(par.at(1),scale*this->alpha,eps);
    EXPECT_NEAR(par.at(2),scale*this->trtAct,eps);
    EXPECT_NEAR(par.at(3),scale*this->trtPre,eps);
  }

  {
    // negative half
    const double scale = -0.5;
    this->m->read();
    this->m->linScale(scale);
    std::vector<double> par = this->m->getPar();
    EXPECT_NEAR(par.at(0),scale*this->intcp,eps);
    EXPECT_NEAR(par.at(1),scale*this->alpha,eps);
    EXPECT_NEAR(par.at(2),scale*this->trtAct,eps);
    EXPECT_NEAR(par.at(3),scale*this->trtPre,eps);
  }
}

TEST_F(TestModel,TestType) {
  EXPECT_EQ(this->m->getType(),INVALID);
  this->m->setType(MLE);
  EXPECT_EQ(this->m->getType(),MLE);
  this->m->setType(MCMC);
  EXPECT_EQ(this->m->getType(),MCMC);
  this->m->setType(MLES);
  EXPECT_EQ(this->m->getType(),MLES);
}

TEST_F(TestModel,TestSetFixSample) {
  EXPECT_EQ(this->m->getFixSample(),0);
  this->m->setFixSample(1);
  EXPECT_EQ(this->m->getFixSample(),1);
}

TEST_F(TestModel,TestEdgeToEdge) {
  EXPECT_FALSE(this->m->getEdgeToEdge());
  this->m->setEdgeToEdge(false);
  EXPECT_FALSE(this->m->getEdgeToEdge());
  this->m->setEdgeToEdge(true);
  EXPECT_TRUE(this->m->getEdgeToEdge());
  this->m->setEdgeToEdge(true);
  EXPECT_TRUE(this->m->getEdgeToEdge());
}

TEST_F(TestModel,TestInfProbs) {
  this->m->read();
  this->m->setFill(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);
  this->m->infProbs(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);

  std::vector<double> infProbs = this->m->infProbs();

  ASSERT_EQ(infProbs.size(),this->system->sD.numNotInfec);

  const double oneNot = 1.0/(1.0 + std::exp(this->intcp - this->alpha*1.0));
  const double twoNot = 1.0/(1.0 + std::exp(this->intcp - this->alpha*2.0));

  EXPECT_NEAR(1.0 - oneNot*twoNot, infProbs.at(0),eps);
  EXPECT_NEAR(1.0 - oneNot*twoNot, infProbs.at(1),eps);

  this->m->setEdgeToEdge(true);
  this->m->setFill(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);
  this->m->infProbs(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);

  infProbs = this->m->infProbs();

  EXPECT_NEAR(1.0 - oneNot, infProbs.at(0),eps);
  EXPECT_NEAR(1.0 - oneNot, infProbs.at(1),eps);
}

TEST_F(TestModel,TestRevProbs) {
  this->m->read();
  this->m->setFill(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);
  this->m->revProbs(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);

  std::vector<double> revProbs = this->m->revProbs();

  ASSERT_EQ(revProbs.size(),this->system->sD.numInfected);

  const double oneNot = 1.0/(1.0 + std::exp(this->intcp - this->alpha*1.0));
  const double twoNot = 1.0/(1.0 + std::exp(this->intcp - this->alpha*2.0));

  EXPECT_NEAR(1.0 - oneNot*twoNot, revProbs.at(0),eps);
  EXPECT_NEAR(1.0 - oneNot*twoNot, revProbs.at(1),eps);

  this->m->setEdgeToEdge(true);
  this->m->setFill(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);
  this->m->revProbs(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);

  revProbs = this->m->revProbs();

  EXPECT_NEAR(1.0 - oneNot, revProbs.at(0),eps);
  EXPECT_NEAR(1.0 - oneNot, revProbs.at(1),eps);
}

TEST_F(TestModel,TestModFill) {
  this->m->read();
  this->m->setFill(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);
  this->m->infProbs(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);

  std::vector<double> infProbs = this->m->infProbs();

  ASSERT_EQ(infProbs.size(),this->system->sD.numNotInfec);

  const double oneNot = 1.0/(1.0 + std::exp(this->intcp - this->alpha*1.0));
  const double twoNot = 1.0/(1.0 + std::exp(this->intcp - this->alpha*2.0));

  EXPECT_NEAR(1.0 - oneNot*twoNot, infProbs.at(0),eps);
  EXPECT_NEAR(1.0 - oneNot*twoNot, infProbs.at(1),eps);

  this->system->tD.a.at(this->system->sD.infected.at(0)) = 1;

  this->m->modFill(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);
  this->m->infProbs(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);

  infProbs = this->m->infProbs();

  const double oneNotTrt = 1.0/(1.0 + std::exp(this->intcp
        - this->alpha*1.0 - this->trtAct));
  const double twoNotTrt = 1.0/(1.0 + std::exp(this->intcp
        - this->alpha*2.0 - this->trtAct));

    EXPECT_NEAR(1.0 - oneNotTrt*twoNot, infProbs.at(0),eps);
    EXPECT_NEAR(1.0 - oneNot*twoNotTrt, infProbs.at(1),eps);
}

TEST_F(TestModel,TestSetQuick) {
  this->m->read();
  this->m->setFill(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);
  this->m->setQuick(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);
  this->m->infProbs(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);

  std::vector<double> infProbs = this->m->infProbs();

  ASSERT_EQ(infProbs.size(),this->system->sD.numNotInfec);

  const double oneNot = 1.0/(1.0 + std::exp(this->intcp - this->alpha*1.0));
  const double twoNot = 1.0/(1.0 + std::exp(this->intcp - this->alpha*2.0));

  EXPECT_NEAR(1.0 - oneNot*twoNot, infProbs.at(0),eps);
  EXPECT_NEAR(1.0 - oneNot*twoNot, infProbs.at(1),eps);

  this->m->setEdgeToEdge(true);
  this->m->setFill(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);
  this->m->setQuick(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);
  this->m->infProbs(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);

  infProbs = this->m->infProbs();

  EXPECT_NEAR(1.0 - oneNot, infProbs.at(0),eps);
  EXPECT_NEAR(1.0 - oneNot, infProbs.at(1),eps);
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

  boost::filesystem::create_directories(njm::sett.srcExt("ParamGDist"));
  njm::toFile(TestModel::intcp,njm::sett.srcExt("ParamGDist/intcp.txt"));
  njm::toFile(TestModel::alpha,njm::sett.srcExt("ParamGDist/alpha.txt"));
  njm::toFile(TestModel::trtAct,njm::sett.srcExt("ParamGDist/trtAct.txt"));
  njm::toFile(TestModel::trtPre,njm::sett.srcExt("ParamGDist/trtPre.txt"));
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
  njm::sett.clean();
  boost::filesystem::remove_all(temp);
  return ret;
}
