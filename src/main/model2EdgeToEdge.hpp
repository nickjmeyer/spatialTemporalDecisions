#ifndef MODEL_2_EDGE_TO_EDGE_HPP__
#define MODEL_2_EDGE_TO_EDGE_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramBeta2.hpp"
#include "paramGravPowGDist.hpp"
#include "paramTrt.hpp"


class Model2EdgeToEdge : public ModelBase {
 protected:
 public:
  Model2EdgeToEdge(){ };
  Model2EdgeToEdge(const FixedData & fD);
  Model2EdgeToEdge(const Model2EdgeToEdge & m);

  virtual Model2EdgeToEdge & operator=(const Model2EdgeToEdge & m);
};


#endif
