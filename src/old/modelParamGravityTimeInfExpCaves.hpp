#ifndef MODEL_PARAM_GRAVITY_TIME_INF_EXP_CAVES_HPP__
#define MODEL_PARAM_GRAVITY_TIME_INF_EXP_CAVES_HPP__


#include <vector>
#include <string>
#include <armadillo>
#include "utilities.hpp"
#include "settings.hpp"
#include "modelParam.hpp"


class GravityTimeInfExpCavesParam : public BaseParam {
 public:
  std::vector<double> infProbs;
  
  std::vector<double> beta;
  double intcp;
  double alpha;
  double power;
  double xi;
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

  virtual std::vector<double> & getInfProbs(){return infProbs;}
  virtual double & getTrtAct(){return trtAct;}
  virtual double & getTrtPre(){return trtPre;}
  virtual arma::mat & getBase(){return infProbsBase;}
  virtual arma::mat & getSep(){return infProbsSep;}
  
};


#endif
