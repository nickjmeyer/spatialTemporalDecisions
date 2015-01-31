#include "modelParamRange.hpp"



void RangeParam::load(){
  njm::fromFile(intcp,njm::sett.srcExt("./RangeParam/intcp.txt"));
  njm::fromFile(range,njm::sett.srcExt("./RangeParam/range.txt"));
  njm::fromFile(alpha,njm::sett.srcExt("./RangeParam/alpha.txt"));
  njm::fromFile(trtAct,njm::sett.srcExt("./RangeParam/trtAct.txt"));
  njm::fromFile(trtPre,njm::sett.srcExt("./RangeParam/trtPre.txt"));
}

std::vector<double> RangeParam::getPar() const {
  std::vector<double> param;
  param.push_back(intcp);
  param.push_back(range);
  param.push_back(alpha);
  param.push_back(trtAct);
  param.push_back(trtPre);
  return param;
}

void RangeParam::putPar(const std::vector<double> & param){
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




inline void RangeParam::setAll(){
  infProbsSep = 1.0/(1.0 + arma::exp(infProbsBase));
}
inline void RangeParam::setRow(const int r){
  infProbsSep.row(r) = 1.0/(1.0 + arma::exp(infProbsBase.row(r)));
}
inline void RangeParam::setCol(const int c){
  infProbsSep.col(c) = 1.0/(1.0 + arma::exp(infProbsBase.col(c)));
}
inline void RangeParam::setInd(const int r, const int c){
  infProbsSep(r,c) = 1.0/(1.0 + std::exp(infProbsBase(r,c)));
}
