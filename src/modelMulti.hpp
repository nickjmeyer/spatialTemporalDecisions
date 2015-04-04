#ifndef MODEL_MULTI_HPP__
#define MODEL_MULTI_HPP__


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include <tuple>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "modelGravityTimeInf.hpp"
#include "modelGravityTimeInfExp.hpp"
#include "modelGravityTimeInfExpCaves.hpp"
#include "modelCave.hpp"
#include "modelRadius.hpp"
#include "modelRange.hpp"


class
MultiModel : public BaseModel {
 public:
  MultiModel();
  MultiModel(const MultiModel & mm);
  ~MultiModel();

  MultiModel & operator=(const MultiModel & mm);
  
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

  virtual BaseParam * getPar(){return m.at(ind)->getPar();}

  virtual void setType(const Estimation & est);

  void fill();

  BaseModel * operator[](const int i);

  int size();
  void modSel(const int & ind);

  
  const static int numModels = 3;

  int ind;

  std::vector<BaseModel *> m;
};



#endif
