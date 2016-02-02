#include "modelParamRadius.hpp"



void RadiusParam::load(){
  njm::fromFile(intcp,njm::sett.srcExt("./RadiusParam/intcp.txt"));
  njm::fromFile(radius,njm::sett.srcExt("./RadiusParam/radius.txt"));
  njm::fromFile(trtAct,njm::sett.srcExt("./RadiusParam/trtAct.txt"));
  njm::fromFile(trtPre,njm::sett.srcExt("./RadiusParam/trtPre.txt"));
}

void RadiusParam::save(){
  njm::toFile(intcp,njm::sett.srcExt("./RadiusParam/intcp.txt"),
	      std::ios_base::out);
  njm::toFile(radius,njm::sett.srcExt("./RadiusParam/radius.txt"),
	      std::ios_base::out);
  njm::toFile(trtAct,njm::sett.srcExt("./RadiusParam/trtAct.txt"),
	      std::ios_base::out);
  njm::toFile(trtPre,njm::sett.srcExt("./RadiusParam/trtPre.txt"),
	      std::ios_base::out);
}

std::vector<double> RadiusParam::getPar() const {
  std::vector<double> param;
  param.push_back(intcp);
  param.push_back(radius);
  param.push_back(trtAct);
  param.push_back(trtPre);
  return param;
}

void RadiusParam::putPar(const std::vector<double> & param){
  std::vector<double>::const_iterator it;
  it = param.end();

  --it;
  trtPre = *it;

  --it;
  trtAct = *it;

  --it;
  radius = *it;

  --it;
  intcp = *it;
}




inline void RadiusParam::setAll(){
  infProbsSep = 1.0/(1.0 + arma::exp(infProbsBase));
}
inline void RadiusParam::setRow(const int r){
  infProbsSep.row(r) = 1.0/(1.0 + arma::exp(infProbsBase.row(r)));
}
inline void RadiusParam::setCol(const int c){
  infProbsSep.col(c) = 1.0/(1.0 + arma::exp(infProbsBase.col(c)));
}
inline void RadiusParam::setInd(const int r, const int c){
  infProbsSep(r,c) = 1.0/(1.0 + std::exp(infProbsBase(r,c)));
}
