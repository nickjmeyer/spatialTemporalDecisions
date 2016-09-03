#ifndef MODEL_GRAVITY_G_DIST_HPP__
#define MODEL_GRAVITY_G_DIST_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramBeta.hpp"
#include "paramGravityGDist.hpp"
#include "paramTrt.hpp"
#include "mcmcGravity.hpp"


class ModelGravityGDist : public ModelBase {
 protected:
 public:
  ModelGravityGDist(){ };
  ModelGravityGDist(const FixedData & fD);
  ModelGravityGDist(const ModelGravityGDist & m);

  virtual ModelGravityGDist & operator=(const ModelGravityGDist & m);

  double tuneTrt(const FixedData & fD);

  GravityMcmc mcmc;
};


#endif
