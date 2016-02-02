#include "modelParamGravityTimeInfSqrt.hpp"


void GravityTimeInfSqrtParam::load(){
  njm::fromFile(intcp,njm::sett.srcExt("./GravityTimeInfSqrtParam/intcp.txt"));
  njm::fromFile(beta,njm::sett.srcExt("./GravityTimeInfSqrtParam/beta.txt"));
  njm::fromFile(alpha,njm::sett.srcExt("./GravityTimeInfSqrtParam/alpha.txt"));
  njm::fromFile(power,njm::sett.srcExt("./GravityTimeInfSqrtParam/power.txt"));
  njm::fromFile(xi,njm::sett.srcExt("./GravityTimeInfSqrtParam/xi.txt"));
  njm::fromFile(trtAct,
		njm::sett.srcExt("./GravityTimeInfSqrtParam/trtAct.txt"));
  njm::fromFile(trtPre,
		njm::sett.srcExt("./GravityTimeInfSqrtParam/trtPre.txt"));
}

void GravityTimeInfSqrtParam::save(){
  njm::toFile(intcp,njm::sett.srcExt("./GravityTimeInfSqrtParam/intcp.txt"),
	      std::ios_base::out);
  njm::toFile(beta,njm::sett.srcExt("./GravityTimeInfSqrtParam/beta.txt"),
	      std::ios_base::out);
  njm::toFile(alpha,njm::sett.srcExt("./GravityTimeInfSqrtParam/alpha.txt"),
	      std::ios_base::out);
  njm::toFile(power,njm::sett.srcExt("./GravityTimeInfSqrtParam/power.txt"),
	      std::ios_base::out);
  njm::toFile(xi,njm::sett.srcExt("./GravityTimeInfSqrtParam/xi.txt"),
	      std::ios_base::out);
  njm::toFile(trtAct,njm::sett.srcExt("./GravityTimeInfSqrtParam/trtAct.txt"),
	      std::ios_base::out);
  njm::toFile(trtPre,njm::sett.srcExt("./GravityTimeInfSqrtParam/trtPre.txt"),
	      std::ios_base::out);
}


std::vector<double> GravityTimeInfSqrtParam::getPar() const {
  std::vector<double> param;
  param.insert(param.end(),beta.begin(),beta.end());
  param.push_back(intcp);
  param.push_back(alpha);
  param.push_back(power);
  param.push_back(xi);
  param.push_back(trtAct);
  param.push_back(trtPre);
  return param;
}

void GravityTimeInfSqrtParam::putPar(const std::vector<double> & param){
  std::vector<double>::const_iterator it;
  it = param.end();

  --it;
  trtPre = *it;

  --it;
  trtAct = *it;

  --it;
  xi = *it;

  --it;
  power = *it;

  --it;
  alpha = *it;

  --it;
  intcp = *it;

  beta.clear();
  beta.insert(beta.begin(),param.begin(),it);
}




inline void GravityTimeInfSqrtParam::setAll(){
  infProbsSep = 1.0/(1.0 + arma::exp(infProbsBase));
}
inline void GravityTimeInfSqrtParam::setRow(const int r){
  infProbsSep.row(r) = 1.0/(1.0 + arma::exp(infProbsBase.row(r)));
}
inline void GravityTimeInfSqrtParam::setCol(const int c){
  infProbsSep.col(c) = 1.0/(1.0 + arma::exp(infProbsBase.col(c)));
}
inline void GravityTimeInfSqrtParam::setInd(const int r, const int c){
  infProbsSep(r,c) = 1.0/(1.0 + std::exp(infProbsBase(r,c)));
}



