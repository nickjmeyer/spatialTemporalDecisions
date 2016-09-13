#ifndef MODEL_INTERCEPT_HPP__
#define MODEL_INTERCEPT_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramTrt.hpp"


class ModelIntercept : public ModelBase {
 protected:
 public:
  ModelIntercept(){ };
  ModelIntercept(const FixedData & fD);
  ModelIntercept(const ModelIntercept & m);

  virtual ModelIntercept & operator=(const ModelIntercept & m);
};


#endif
