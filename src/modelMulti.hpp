#ifndef MODEL_MULTI_HPP__
#define MODEL_MULTI_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include <tuple>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "modelCave.hpp"
#include "modelRadius.hpp"
#include "modelRange.hpp"


class
MultiModel : public BaseModel {
 public:
  virtual void load(const SimData & sD,
		    const TrtData & tD,
		    const FixedData & fD,
		    const DynamicData & dD);

  virtual void infProbs(const SimData & sD,
			const TrtData & tD,
			const FixedData & fD,
			const DynamicData & dD);
  
  virtual void update(const SimData & sD,
		      const TrtData & tD,
		      const FixedData & fD,
		      const DynamicData & dD);
  
  virtual double oneOnOne(const int notNode, const int infNode,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD) const;

  virtual void fit(const SimData & sD, const TrtData & tD,
		   const FixedData & fD, const DynamicData & dD,
		   const int & useInit);
  
  void modSel(const int & ind);
  int ind;

  std::vector<BaseModel> m;
};



#endif
