#include "modelParamGravityTimeInfLog.hpp"


void GravityTimeInfLogParam::load(){
  njm::fromFile(intcp,njm::sett.srcExt("./GravityTimeInfLogParam/intcp.txt"));
  njm::fromFile(beta,njm::sett.srcExt("./GravityTimeInfLogParam/beta.txt"));
  njm::fromFile(alpha,njm::sett.srcExt("./GravityTimeInfLogParam/alpha.txt"));
  njm::fromFile(power,njm::sett.srcExt("./GravityTimeInfLogParam/power.txt"));
  njm::fromFile(xi,njm::sett.srcExt("./GravityTimeInfLogParam/xi.txt"));
  njm::fromFile(trtAct,
		njm::sett.srcExt("./GravityTimeInfLogParam/trtAct.txt"));
  njm::fromFile(trtPre,
		njm::sett.srcExt("./GravityTimeInfLogParam/trtPre.txt"));
}

void GravityTimeInfLogParam::save(){
  njm::toFile(intcp,njm::sett.srcExt("./GravityTimeInfLogParam/intcp.txt"),
	      std::ios_base::out);
  njm::toFile(beta,njm::sett.srcExt("./GravityTimeInfLogParam/beta.txt"),
	      std::ios_base::out);
  njm::toFile(alpha,njm::sett.srcExt("./GravityTimeInfLogParam/alpha.txt"),
	      std::ios_base::out);
  njm::toFile(power,njm::sett.srcExt("./GravityTimeInfLogParam/power.txt"),
	      std::ios_base::out);
  njm::toFile(xi,njm::sett.srcExt("./GravityTimeInfLogParam/xi.txt"),
	      std::ios_base::out);
  njm::toFile(trtAct,njm::sett.srcExt("./GravityTimeInfLogParam/trtAct.txt"),
	      std::ios_base::out);
  njm::toFile(trtPre,njm::sett.srcExt("./GravityTimeInfLogParam/trtPre.txt"),
	      std::ios_base::out);
}


std::vector<double> GravityTimeInfLogParam::getPar() const {
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

void GravityTimeInfLogParam::putPar(const std::vector<double> & param){
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




inline void GravityTimeInfLogParam::setAll(){
  infProbsSep = 1.0/(1.0 + arma::exp(infProbsBase));
}
inline void GravityTimeInfLogParam::setRow(const int r){
  infProbsSep.row(r) = 1.0/(1.0 + arma::exp(infProbsBase.row(r)));
}
inline void GravityTimeInfLogParam::setCol(const int c){
  infProbsSep.col(c) = 1.0/(1.0 + arma::exp(infProbsBase.col(c)));
}
inline void GravityTimeInfLogParam::setInd(const int r, const int c){
  infProbsSep(r,c) = 1.0/(1.0 + std::exp(infProbsBase(r,c)));
}



