#ifndef MODEL_E_DIST_HPP
#define MODEL_E_DIST_HPP


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramEDist.hpp"
#include "paramTrt.hpp"


class ModelEDist : public ModelBase {
 protected:
 public:
  ModelEDist();
  ModelEDist(const FixedData & fD);
  ModelEDist(const ModelEDist & m);

  virtual ModelEDist & operator=(const ModelEDist & m);
};


#endif
