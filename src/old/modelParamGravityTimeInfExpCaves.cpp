#include "modelParamGravityTimeInfExpCaves.hpp"


void GravityTimeInfExpCavesParam::load(){
  njm::fromFile(intcp,
		njm::sett.srcExt("./GravityTimeInfExpCavesParam/intcp.txt"));
  njm::fromFile(beta,
		njm::sett.srcExt("./GravityTimeInfExpCavesParam/beta.txt"));
  njm::fromFile(alpha,
		njm::sett.srcExt("./GravityTimeInfExpCavesParam/alpha.txt"));
  njm::fromFile(power,
		njm::sett.srcExt("./GravityTimeInfExpCavesParam/power.txt"));
  njm::fromFile(xi,
		njm::sett.srcExt("./GravityTimeInfExpCavesParam/xi.txt"));
  njm::fromFile(trtAct,
		njm::sett.srcExt("./GravityTimeInfExpCavesParam/trtAct.txt"));
  njm::fromFile(trtPre,
		njm::sett.srcExt("./GravityTimeInfExpCavesParam/trtPre.txt"));
}

void GravityTimeInfExpCavesParam::save(){
  njm::toFile(intcp,
	      njm::sett.srcExt("./GravityTimeInfExpCavesParam/intcp.txt"),
	      std::ios_base::out);
  njm::toFile(beta,
	      njm::sett.srcExt("./GravityTimeInfExpCavesParam/beta.txt"),
	      std::ios_base::out);
  njm::toFile(alpha,
	      njm::sett.srcExt("./GravityTimeInfExpCavesParam/alpha.txt"),
	      std::ios_base::out);
  njm::toFile(power,
	      njm::sett.srcExt("./GravityTimeInfExpCavesParam/power.txt"),
	      std::ios_base::out);
  njm::toFile(xi,
	      njm::sett.srcExt("./GravityTimeInfExpCavesParam/xi.txt"),
	      std::ios_base::out);
  njm::toFile(trtAct,
	      njm::sett.srcExt("./GravityTimeInfExpCavesParam/trtAct.txt"),
	      std::ios_base::out);
  njm::toFile(trtPre,
	      njm::sett.srcExt("./GravityTimeInfExpCavesParam/trtPre.txt"),
	      std::ios_base::out);
}


std::vector<double> GravityTimeInfExpCavesParam::getPar() const {
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

void GravityTimeInfExpCavesParam::putPar(const std::vector<double> & param){
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




inline void GravityTimeInfExpCavesParam::setAll(){
  infProbsSep = 1.0/(1.0 + arma::exp(infProbsBase));
}
inline void GravityTimeInfExpCavesParam::setRow(const int r){
  infProbsSep.row(r) = 1.0/(1.0 + arma::exp(infProbsBase.row(r)));
}
inline void GravityTimeInfExpCavesParam::setCol(const int c){
  infProbsSep.col(c) = 1.0/(1.0 + arma::exp(infProbsBase.col(c)));
}
inline void GravityTimeInfExpCavesParam::setInd(const int r, const int c){
  infProbsSep(r,c) = 1.0/(1.0 + std::exp(infProbsBase(r,c)));
}



