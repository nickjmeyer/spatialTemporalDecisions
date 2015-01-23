#include "modelParam.hpp"


void GravityParam::load(){
  njm::fromFile(intcp,njm::sett.srcExt("./GravityParam/intcp.txt"));
  njm::fromFile(beta,njm::sett.srcExt("./GravityParam/beta.txt"));
  njm::fromFile(alpha,njm::sett.srcExt("./GravityParam/alpha.txt"));
  njm::fromFile(power,njm::sett.srcExt("./GravityParam/power.txt"));
  njm::fromFile(trtAct,njm::sett.srcExt("./GravityParam/trtAct.txt"));
  njm::fromFile(trtPre,njm::sett.srcExt("./GravityParam/trtPre.txt"));
}

std::vector<double> GravityParam::getPar() const {
  std::vector<double> param;
  param.insert(param.end(),beta.begin(),beta.end());
  param.push_back(intcp);
  param.push_back(alpha);
  param.push_back(power);
  param.push_back(trtAct);
  param.push_back(trtPre);
  return param;
}

void GravityParam::putPar(const std::vector<double> & param){
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

  beta.clear();
  beta.insert(beta.begin(),param.begin(),it);
}




inline void GravityParam::setAll(){
  infProbsSep = 1.0/(1.0 + arma::exp(infProbsBase));
}
inline void GravityParam::setRow(const int r){
  infProbsSep.row(r) = 1.0/(1.0 + arma::exp(infProbsBase.row(r)));
}
inline void GravityParam::setCol(const int c){
  infProbsSep.col(c) = 1.0/(1.0 + arma::exp(infProbsBase.col(c)));
}
inline void GravityParam::setInd(const int r, const int c){
  infProbsSep(r,c) = 1.0/(1.0 + std::exp(infProbsBase(r,c)));
}



