#include "modelParamCave.hpp"



void CaveParam::load(){
  njm::fromFile(intcp,njm::sett.srcExt("./CaveParam/intcp.txt"));
  njm::fromFile(cave,njm::sett.srcExt("./CaveParam/cave.txt"));
  njm::fromFile(trtAct,njm::sett.srcExt("./CaveParam/trtAct.txt"));
  njm::fromFile(trtPre,njm::sett.srcExt("./CaveParam/trtPre.txt"));
}

std::vector<double> CaveParam::getPar() const {
  std::vector<double> param;
  param.push_back(intcp);
  param.push_back(cave);
  param.push_back(trtAct);
  param.push_back(trtPre);
  return param;
}

void CaveParam::putPar(const std::vector<double> & param){
  std::vector<double>::const_iterator it;
  it = param.end();

  --it;
  trtPre = *it;

  --it;
  trtAct = *it;

  --it;
  cave = *it;

  --it;
  intcp = *it;
}




inline void CaveParam::setAll(){
  infProbsSep = 1.0/(1.0 + arma::exp(infProbsBase));
}
inline void CaveParam::setRow(const int r){
  infProbsSep.row(r) = 1.0/(1.0 + arma::exp(infProbsBase.row(r)));
}
inline void CaveParam::setCol(const int c){
  infProbsSep.col(c) = 1.0/(1.0 + arma::exp(infProbsBase.col(c)));
}
inline void CaveParam::setInd(const int r, const int c){
  infProbsSep(r,c) = 1.0/(1.0 + std::exp(infProbsBase(r,c)));
}
