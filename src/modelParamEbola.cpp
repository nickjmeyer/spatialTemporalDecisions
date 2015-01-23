#include "modelParamEbola.hpp"



void EbolaParam::load(){
  njm::fromFile(intcp,njm::sett.srcExt("./EbolaParam/intcp.txt"));
  njm::fromFile(alpha,njm::sett.srcExt("./EbolaParam/alpha.txt"));
  njm::fromFile(power,njm::sett.srcExt("./EbolaParam/power.txt"));
  njm::fromFile(trtAct,njm::sett.srcExt("./EbolaParam/trtAct.txt"));
  njm::fromFile(trtPre,njm::sett.srcExt("./EbolaParam/trtPre.txt"));
}

std::vector<double> EbolaParam::getPar() const {
  std::vector<double> param;
  param.push_back(intcp);
  param.push_back(alpha);
  param.push_back(power);
  param.push_back(trtAct);
  param.push_back(trtPre);
  return param;
}

void EbolaParam::putPar(const std::vector<double> & param){
  std::vector<double>::const_iterator it;
  it = param.end();

  --it;
  trtPre = *it;

  --it;
  trtAct = *it;

  --it;
  power = *it;

  --it;
  alpha = *it;

  --it;
  intcp = *it;
}




inline void EbolaParam::setAll(){
  infProbsSep = 1.0/(1.0 + arma::exp(infProbsBase));
}
inline void EbolaParam::setRow(const int r){
  infProbsSep.row(r) = 1.0/(1.0 + arma::exp(infProbsBase.row(r)));
}
inline void EbolaParam::setCol(const int c){
  infProbsSep.col(c) = 1.0/(1.0 + arma::exp(infProbsBase.col(c)));
}
inline void EbolaParam::setInd(const int r, const int c){
  infProbsSep(r,c) = 1.0/(1.0 + std::exp(infProbsBase(r,c)));
}
