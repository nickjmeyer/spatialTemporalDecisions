#ifndef MODEL_PARAM_CAVE_HPP__
#define MODEL_PARAM_CAVE_HPP__


#include <vector>
#include <string>
#include <armadillo>
#include "utilities.hpp"
#include "modelParam.hpp"
#include "settings.hpp"

class CaveParam : public BaseParam {
 public:
  std::vector<double> infProbs;
  
  double intcp;
  double cave;
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

  virtual double & getTrtAct(){return trtAct;}
  virtual double & getTrtPre(){return trtPre;}
  virtual arma::mat & getBase(){return infProbsBase;}
  virtual arma::mat & getSep(){return infProbsSep;}
  
};


#endif
