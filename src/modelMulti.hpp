#ifndef MODEL_MULTI_HPP__
#define MODEL_MULTI_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include <tuple>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "modelParamMulti.hpp"
#include "modelCave.hpp"
#include "modelRadius.hpp"
#include "modelRange.hpp"


class
MultiModel : public BaseModel {
 public:
  virtual double oneOnOne(const int notNode, const int infNode,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD) const;

  
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD);
  void fit(const SimData & sD, const TrtData & tD,
	   const FixedData & fD, const DynamicData & dD,
	   const std::vector<double> & mPV);

  void modSel(const int & ind);
  int ind;

  std::vector<BaseModel> m;

  Estimation fitType;
};



#endif
