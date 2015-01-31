#ifndef MODEL_PARAM_RANGE_HPP__
#define MODEL_PARAM_RANGE_HPP__


#include <vector>
#include <string>
#include <armadillo>
#include "utilities.hpp"
#include "modelParam.hpp"
#include "settings.hpp"

class RangeParam : public BaseParam {
 public:
  std::vector<double> infProbs;
  
  std::vector<double> beta;
  double intcp;
  double range;
  double alpha;
  double trtAct;
  double trtPre;

  arma::mat infProbsBase;
  arma::mat infProbsSep;

  virtual void load();
  
  virtual std::vector<double> getPar() const;

  virtual void putPar(const std::vector<double> & param);

  virtual void setAll();
  virtual void setRow(const int r);
  virtual void setCol(const int c);
  virtual void setInd(const int r, const int c);
};


#endif
