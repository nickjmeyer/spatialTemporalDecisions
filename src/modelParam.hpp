#ifndef MODEL_PARAM_HPP__
#define MODEL_PARAM_HPP__


#include <vector>
#include <string>
#include <armadillo>
#include "utilities.hpp"
#include "settings.hpp"


class BaseParam {
 public:
  std::vector<double> infProbs;
  
  virtual void load() = 0;
  virtual std::vector<double> getPar() const = 0;
  virtual void putPar(const std::vector<double> & param) = 0;
};


class GravityParam : public BaseParam {
 public:
  std::vector<double> infProbs;
  
  std::vector<double> beta;
  double intcp;
  double alpha;
  double power;
  double trtAct;
  double trtPre;

  arma::mat infProbsBase;
  arma::mat infProbsSep;

  virtual void load();
  virtual void save();
  
  virtual std::vector<double> getPar() const;

  virtual void putPar(const std::vector<double> & param);

  virtual void setAll();
  virtual void setRow(const int r);
  virtual void setCol(const int c);
  virtual void setInd(const int r, const int c);
};


#endif
