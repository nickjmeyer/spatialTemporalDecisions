#ifndef MODEL_TIME_EXP_CAVES_G_POW_G_DIST_TREND_POW_CON_HPP__
#define MODEL_TIME_EXP_CAVES_G_POW_G_DIST_TREND_POW_CON_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "timer.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramBeta.hpp"
#include "paramGravPowGDist.hpp"
#include "paramTrendPowCon.hpp"
#include "paramTimeExpCaves.hpp"
#include "paramTrt.hpp"
#include "mcmcGravityTimeInfExpCavesTrendPowCon.hpp"


class ModelTimeExpCavesGPowGDistTrendPowCon : public ModelBase {
 protected:
 public:
  ModelTimeExpCavesGPowGDistTrendPowCon(){ };
  ModelTimeExpCavesGPowGDistTrendPowCon(const FixedData & fD);
  ModelTimeExpCavesGPowGDistTrendPowCon
  (const ModelTimeExpCavesGPowGDistTrendPowCon & m);

  virtual ModelTimeExpCavesGPowGDistTrendPowCon &
  operator=(const ModelTimeExpCavesGPowGDistTrendPowCon & m);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   std::vector<double> all);
};


class ModelTimeExpCavesGPowGDistTrendPowConFitData {
 public:
  ModelTimeExpCavesGPowGDistTrendPowConFitData
  (const ModelTimeExpCavesGPowGDistTrendPowCon & m,
   const std::vector<double> & all,
   const FixedData & fD,
   const
   std::vector<std::vector<int> > & history);

  ModelTimeExpCavesGPowGDistTrendPowCon m;
  FixedData fD;
  std::vector<std::vector<int> > history;
  std::vector<std::vector<int> > timeInf;
};


double
modelTimeExpCavesGPowGDistTrendPowConFitObjFn (const gsl_vector * x,
					       void * params);




#endif
