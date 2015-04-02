#ifndef MODEL_PARAM_MULTI_HPP__
#define MODEL_PARAM_MULTI_HPP__


#include <vector>
#include <string>
#include <armadillo>
#include "utilities.hpp"
#include "settings.hpp"
#include "modelParam.hpp"
#include "modelParamCave.hpp"
#include "modelParamRadius.hpp"
#include "modelParamRange.hpp"


class MultiParam : public BaseParam {
 public:

  virtual void load();
  virtual void save();
  
  virtual std::vector<double> getPar() const;

  virtual void putPar(const std::vector<double> & param);

  virtual void setAll();
  virtual void setRow(const int r);
  virtual void setCol(const int c);
  virtual void setInd(const int r, const int c);
  
  std::tuple<RadiusParam,RangeParam,CaveParam> pars;
};


#endif
