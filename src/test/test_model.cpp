#include <gflags/gflags.h>
#include <gtest/gtest.h>
#include <glog/logging.h>
#include <boost/filesystem.hpp>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include "settings.hpp"
#include "data.hpp"
#include "system.hpp"
#include "paramIntercept.hpp"

const float eps = 1e-6;

const std::vector<std::string> modNames = {
  "2GPowGDist",
  "2GravityGDist",
  "2GravityEDist",
  "GDist",
  "EDist",
  "GravityGDist",
  "Intercept"};

const std::vector<std::string> parNames = {
  "intcp",
  "beta0",
  "beta1",
  "beta2",
  "beta3",
  "alpha",
  "power",
  "gPow",
  "trtAct",
  "trtPre"};

template <typename T>
class GradientChecker {
public:
  System<T,T> * system;
  T * m;
  int gradientVar;
};

template <typename T>
class HessianChecker {
public:
  System<T,T> * system;
  T * m;
  int gradientVar;
  int hessianVar;
};

template <typename T>
double f (double x, void * params) {
  GradientChecker<T> * gc = static_cast<GradientChecker<T>*>(params);
  std::vector<double> par = gc->m->getPar();
  par.at(gc->gradientVar) = x;
  gc->m->putPar(par.begin());
  return gc->m->logll(gc->system->sD,gc->system->tD,
    gc->system->fD,gc->system->dD);
}

template <typename T>
double fGrad (double x, void * params) {
  HessianChecker<T> * hc = static_cast<HessianChecker<T>*>(params);
  std::vector<double> par = hc->m->getPar();
  par.at(hc->hessianVar) = x;
  hc->m->putPar(par.begin());
  return hc->m->logllGrad(hc->system->sD,hc->system->tD,
    hc->system->fD,hc->system->dD).at(hc->gradientVar);
}



template <typename T>
class TestModel : public ::testing::Test {
public:
  T * m;

  System<T,T> * system;

  TestModel() {
    system = new System<T,T>();
    system->reset({0});
    std::vector<double> mockProbs = {1.0,0.0,0.0};
    system->nextPoint(mockProbs);

    m = new T(system->fD);
  }

  ~TestModel() {
    delete m;
    delete system;
  }
};


typedef ::testing::Types<
    ModelIntercept,
    ModelGDist,
    ModelEDist,
    ModelGravityGDist,
    Model2GPowGDist,
    Model2GravityGDist,
    Model2GravityEDist
    > MyTypes;
TYPED_TEST_CASE(TestModel,MyTypes);

TYPED_TEST(TestModel,TestInit) {
  for (int i = 0; i < this->m->numPars; ++i) {
    EXPECT_NEAR(this->m->getPar().at(i),0.0,eps);
  }

  unsigned int offset = 0;
  for (int i = 0; i < this->m->pars.size(); ++i) {
      EXPECT_EQ(offset,this->m->pars[i]->getOffset());
      EXPECT_EQ(this->m->numPars,this->m->pars[i]->getTotNumPars());
      offset += this->m->pars[i]->size();
  }
  EXPECT_EQ(offset,this->m->numPars);
}

TYPED_TEST(TestModel,TestReadSave) {
  this->m->read();

  std::vector<double> oldPar = this->m->getPar();
  std::vector<double> newPar(this->m->numPars,0);
  for (int i = 0; i < this->m->numPars; ++i) {
    newPar.at(i) = njm::runif(-0.5,0.5);
  }
  this->m->putPar(newPar.begin());
  this->m->save();

  std::vector<double> zeroPar(this->m->numPars,0);
  this->m->putPar(zeroPar.begin());
  for (int i = 0; i < this->m->numPars; ++i) {
    EXPECT_EQ(this->m->getPar().at(i),0.);
  }

  this->m->read();
  for (int i = 0; i < this->m->numPars; ++i) {
    EXPECT_NEAR(this->m->getPar().at(i),newPar.at(i),eps);
  }

  // now make sure the old parameters write properly because other
  // test cases might depend on it
  this->m->putPar(oldPar.begin());
  this->m->save();
  this->m->putPar(zeroPar.begin());
  for (int i = 0; i < this->m->numPars; ++i) {
    EXPECT_EQ(this->m->getPar().at(i),0.);
  }
  this->m->read();
  for (int i = 0; i < this->m->numPars; ++i) {
    EXPECT_NEAR(this->m->getPar().at(i),oldPar.at(i),eps);
  }
}

TYPED_TEST(TestModel,TestLinScale) {
  this->m->setFill(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);
  const std::vector<double> probs = this->m->probs;
  EXPECT_EQ(probs.size(),this->system->fD.numNodes*this->system->fD.numNodes);

  const std::vector<double> scaleVals = {1.0,0.5,2.0,
                                         -1.0,-0.5,-2.0};
  for (int i = 0; i < scaleVals.size(); ++i) {
    const double scale = scaleVals.at(i);
    this->m->linScale(scale);

    const std::vector<double> scaleProbs = this->m->probs;
    EXPECT_EQ(scaleProbs.size(),probs.size());
    for (int j = 0; j < scaleProbs.size(); ++j) {
      EXPECT_NEAR(scaleProbs.at(j),probs.at(j)*scale,eps);
    }
  }
}

TYPED_TEST(TestModel,TestType) {
  EXPECT_EQ(this->m->getType(),INVALID);
  this->m->setType(MLE);
  EXPECT_EQ(this->m->getType(),MLE);
  this->m->setType(MCMC);
  EXPECT_EQ(this->m->getType(),MCMC);
  this->m->setType(MLES);
  EXPECT_EQ(this->m->getType(),MLES);
}

TYPED_TEST(TestModel,TestSetFixSample) {
  EXPECT_EQ(this->m->getFixSample(),0);
  this->m->setFixSample(1);
  EXPECT_EQ(this->m->getFixSample(),1);
}

TYPED_TEST(TestModel,TestEdgeToEdge) {
  EXPECT_FALSE(this->m->getEdgeToEdge());
  this->m->setEdgeToEdge(false);
  EXPECT_FALSE(this->m->getEdgeToEdge());
  this->m->setEdgeToEdge(true);
  EXPECT_TRUE(this->m->getEdgeToEdge());
  this->m->setEdgeToEdge(true);
  EXPECT_TRUE(this->m->getEdgeToEdge());
}

TYPED_TEST(TestModel,TestInfProbs) {
  this->m->read();
  this->m->setFill(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);
  this->m->infProbs(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);

  std::vector<double> infProbs = this->m->infProbs();

  ASSERT_EQ(infProbs.size(),this->system->sD.numNotInfec);
  ASSERT_EQ(this->system->sD.numNotInfec,2);
  ASSERT_EQ(this->system->sD.infected.size(),2);
  ASSERT_EQ(this->system->sD.numInfected,2);
  ASSERT_EQ(this->system->sD.notInfec.size(),2);
  ASSERT_EQ(this->system->sD.infected.at(0),0);
  ASSERT_EQ(this->system->sD.infected.at(1),1);
  ASSERT_EQ(this->system->sD.notInfec.at(0),2);
  ASSERT_EQ(this->system->sD.notInfec.at(1),3);

  const double prob0To2 = 1.0/(1.0 +
    std::exp(this->m->probs[2*this->system->fD.numNodes + 0]));
  const double prob0To3 = 1.0/(1.0 +
    std::exp(this->m->probs[3*this->system->fD.numNodes + 0]));

  const double prob1To2 = 1.0/(1.0 +
    std::exp(this->m->probs[2*this->system->fD.numNodes + 1]));
  const double prob1To3 = 1.0/(1.0 +
    std::exp(this->m->probs[3*this->system->fD.numNodes + 1]));

  EXPECT_NEAR(1.0 - prob0To2*prob1To2, infProbs.at(0),eps);
  EXPECT_NEAR(1.0 - prob0To3*prob1To3, infProbs.at(1),eps);

  this->m->setEdgeToEdge(true);
  this->m->setFill(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);
  this->m->infProbs(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);

  infProbs = this->m->infProbs();

  EXPECT_NEAR(1.0 - prob0To2, infProbs.at(0),eps);
  EXPECT_NEAR(1.0 - prob1To3, infProbs.at(1),eps);
}

TYPED_TEST(TestModel,TestRevProbs) {
  this->m->read();
  this->m->setFill(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);
  this->m->revProbs(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);

  std::vector<double> revProbs = this->m->revProbs();

  ASSERT_EQ(revProbs.size(),this->system->sD.numNotInfec);
  ASSERT_EQ(this->system->sD.numNotInfec,2);
  ASSERT_EQ(this->system->sD.infected.size(),2);
  ASSERT_EQ(this->system->sD.numInfected,2);
  ASSERT_EQ(this->system->sD.notInfec.size(),2);
  ASSERT_EQ(this->system->sD.infected.at(0),0);
  ASSERT_EQ(this->system->sD.infected.at(1),1);
  ASSERT_EQ(this->system->sD.notInfec.at(0),2);
  ASSERT_EQ(this->system->sD.notInfec.at(1),3);

  const double prob2To0 = 1.0/(1.0 +
    std::exp(this->m->probs[0*this->system->fD.numNodes + 2]));
  const double prob3To0 = 1.0/(1.0 +
    std::exp(this->m->probs[0*this->system->fD.numNodes + 3]));

  const double prob2To1 = 1.0/(1.0 +
    std::exp(this->m->probs[1*this->system->fD.numNodes + 2]));
  const double prob3To1 = 1.0/(1.0 +
    std::exp(this->m->probs[1*this->system->fD.numNodes + 3]));

  EXPECT_NEAR(1.0 - prob2To0*prob3To0, revProbs.at(0),eps);
  EXPECT_NEAR(1.0 - prob2To1*prob3To1, revProbs.at(1),eps);

  this->m->setEdgeToEdge(true);
  this->m->setFill(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);
  this->m->revProbs(this->system->sD,this->system->tD,
    this->system->fD,this->system->dD);

  revProbs = this->m->revProbs();

  EXPECT_NEAR(1.0 - prob2To0, revProbs.at(0),eps);
  EXPECT_NEAR(1.0 - prob3To1, revProbs.at(1),eps);
}

TYPED_TEST(TestModel,TestModFill) {
  // this->m->read();
  // this->m->setFill(this->system->sD,this->system->tD,
  //   this->system->fD,this->system->dD);
  // this->m->infProbs(this->system->sD,this->system->tD,
  //   this->system->fD,this->system->dD);

  // std::vector<double> infProbs = this->m->infProbs();

  // ASSERT_EQ(infProbs.size(),this->system->sD.numNotInfec);

  // const double oneNot = 1.0/(1.0 + std::exp(this->intcp - this->alpha*1.0));
  // const double twoNot = 1.0/(1.0 + std::exp(this->intcp - this->alpha*2.0));

  // EXPECT_NEAR(1.0 - oneNot*twoNot, infProbs.at(0),eps);
  // EXPECT_NEAR(1.0 - oneNot*twoNot, infProbs.at(1),eps);

  // this->system->tD.a.at(this->system->sD.infected.at(0)) = 1;

  // this->m->modFill(this->system->sD,this->system->tD,
  //   this->system->fD,this->system->dD);
  // this->m->infProbs(this->system->sD,this->system->tD,
  //   this->system->fD,this->system->dD);

  // infProbs = this->m->infProbs();

  // const double oneNotTrt = 1.0/(1.0 + std::exp(this->intcp
  //       - this->alpha*1.0 - this->trtAct));
  // const double twoNotTrt = 1.0/(1.0 + std::exp(this->intcp
  //       - this->alpha*2.0 - this->trtAct));

  //   EXPECT_NEAR(1.0 - oneNotTrt*twoNot, infProbs.at(0),eps);
  //   EXPECT_NEAR(1.0 - oneNot*twoNotTrt, infProbs.at(1),eps);
}

TYPED_TEST(TestModel,TestSetQuick) {
  // this->m->read();
  // this->m->setFill(this->system->sD,this->system->tD,
  //   this->system->fD,this->system->dD);
  // this->m->setQuick(this->system->sD,this->system->tD,
  //   this->system->fD,this->system->dD);
  // this->m->infProbs(this->system->sD,this->system->tD,
  //   this->system->fD,this->system->dD);

  // std::vector<double> infProbs = this->m->infProbs();

  // ASSERT_EQ(infProbs.size(),this->system->sD.numNotInfec);

  // const double oneNot = 1.0/(1.0 + std::exp(this->intcp - this->alpha*1.0));
  // const double twoNot = 1.0/(1.0 + std::exp(this->intcp - this->alpha*2.0));

  // EXPECT_NEAR(1.0 - oneNot*twoNot, infProbs.at(0),eps);
  // EXPECT_NEAR(1.0 - oneNot*twoNot, infProbs.at(1),eps);

  // this->m->setEdgeToEdge(true);
  // this->m->setFill(this->system->sD,this->system->tD,
  //   this->system->fD,this->system->dD);
  // this->m->setQuick(this->system->sD,this->system->tD,
  //   this->system->fD,this->system->dD);
  // this->m->infProbs(this->system->sD,this->system->tD,
  //   this->system->fD,this->system->dD);

  // infProbs = this->m->infProbs();

  // EXPECT_NEAR(1.0 - oneNot, infProbs.at(0),eps);
  // EXPECT_NEAR(1.0 - oneNot, infProbs.at(1),eps);
}

TYPED_TEST(TestModel,TestSetPar) {
}

TYPED_TEST(TestModel,TestPcPartial) {
    this->m->read();
    this->m->setFill(this->system->sD,this->system->tD,
            this->system->fD,this->system->dD,true);
    for (int i = 0; i < this->system->sD.numNotInfec; i++) {
        const int nNode = this->system->sD.notInfec.at(i);

        for (int j = 0; j < this->system->sD.numInfected; ++j) {
            const int iNode = this->system->sD.infected.at(j);

            for (int k = 0; k < this->m->numPars; ++k) {
                const int pcPartialIndex =
                    nNode * this->system->fD.numNodes * this->m->numPars +
                    iNode * this->m->numPars +
                    k;
                EXPECT_EQ(this->m->pcPartial.at(pcPartialIndex),
                        this->m->partial(nNode,iNode,this->system->sD,
                                this->system->tD,this->system->fD,
                                this->system->dD).at(k))
                    << "failed for not infected node " << nNode
                    << " and infected node " << iNode
                    << " for param " << k;
            }
        }
    }

    std::vector<double> infProbs(this->system->sD.numNotInfec,0.0);
    infProbs.at(0) = 1.0;
    this->system->nextPoint(infProbs);

    this->system->tD.a.at(this->system->sD.infected.at(0)) = 1;
    this->system->tD.p.at(this->system->sD.notInfec.at(0)) = 1;

    this->m->modFill(this->system->sD,this->system->tD,
            this->system->fD,this->system->dD,true);
    for (int i = 0; i < this->system->sD.numNotInfec; i++) {
        const int nNode = this->system->sD.notInfec.at(i);

        for (int j = 0; j < this->system->sD.numInfected; ++j) {
            const int iNode = this->system->sD.infected.at(j);

            for (int k = 0; k < this->m->numPars; ++k) {
                const int pcPartialIndex =
                    nNode * this->system->fD.numNodes * this->m->numPars +
                    iNode * this->m->numPars +
                    k;
                EXPECT_EQ(this->m->pcPartial.at(pcPartialIndex),
                        this->m->partial(nNode,iNode,this->system->sD,
                                this->system->tD,this->system->fD,
                                this->system->dD).at(k))
                    << "failed for not infected node " << nNode
                    << " and infected node " << iNode
                    << " for param " << k;
            }
        }
    }
}


TYPED_TEST(TestModel,TestLogllGrad) {
  this->m->read();

  std::vector<double> par = this->m->getPar();

  const double val = this->m->logll(this->system->sD,this->system->tD,
          this->system->fD,this->system->dD);

  const std::vector<double> gradVal = this->m->logllGrad(
          this->system->sD,this->system->tD,
          this->system->fD,this->system->dD);

  const std::pair<double,std::vector<double> > both = this->m->logllBoth(
          this->system->sD,this->system->tD,
          this->system->fD,this->system->dD);

  EXPECT_EQ(val,both.first);

  for (int i = 0; i < this->m->numPars; ++i) {
    this->m->putPar(par.begin());

    GradientChecker<TypeParam> gc;
    gc.system = this->system;
    gc.m = this->m;
    gc.gradientVar = i;

    gsl_function F;
    F.function = &f<TypeParam>;
    F.params = &gc;

    double result;
    double abserr;
    gsl_deriv_central(&F,par.at(i),1e-8,&result,&abserr);

    EXPECT_EQ(gradVal.at(i),both.second.at(i));
    EXPECT_NEAR(result,gradVal.at(i),eps);
    EXPECT_NEAR(result,both.second.at(i),eps);
  }
}

TYPED_TEST(TestModel,TestLogllHess) {
  this->m->read();


  std::vector<double> par = this->m->getPar();

  for (int i = 0; i < this->m->numPars; ++i) {
    for (int j = 0; j < this->m->numPars; ++j) {
      this->m->putPar(par.begin());

      const double val0 = this->m->logllHess(this->system->sD,this->system->tD,
        this->system->fD,this->system->dD).at(i*this->m->numPars + j);
      const double val1 = this->m->logllHess(this->system->sD,this->system->tD,
        this->system->fD,this->system->dD).at(j*this->m->numPars + i);
      EXPECT_NEAR(val0,val1,eps);

      HessianChecker<TypeParam> hc;
      hc.system = this->system;
      hc.m = this->m;
      hc.gradientVar = i;
      hc.hessianVar = j;

      gsl_function F;
      F.function = &fGrad<TypeParam>;
      F.params = &hc;

      double result;
      double abserr;
      gsl_deriv_central(&F,par.at(j),1e-8,&result,&abserr);

      EXPECT_NEAR(result,val0,eps)
        << "Gradient var is " << i << " and Hessian var is " << j;
    }
  }
}

TYPED_TEST(TestModel,TestSetFisher) {
}

TYPED_TEST(TestModel,TestSample) {
}

TYPED_TEST(TestModel,TestRevert) {
}

TYPED_TEST(TestModel,TestFit) {
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

  for(int i = 0; i < modNames.size(); ++i) {
    const std::string modName = modNames.at(i);
    boost::filesystem::create_directories(njm::sett.srcExt("Param"+modName));
    // some pars aren't used in every model, but its cheap just to
    // write all parameters for all models rather than listing them
    // separately for each model
    for (int j = 0; j < parNames.size(); ++j) {
      const std::string parName = parNames.at(j);
      njm::toFile(njm::runif(-0.5,0.5),njm::sett.srcExt("Param"+modName+
          "/"+parName+".txt"));
    }

  }
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
