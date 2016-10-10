#include "runStats.hpp"

RunStats::RunStats(){
  n = 0U;
  vals.clear();
  sMean_ = 0.0;
  sSqMean_ = 0.0;
  sVar_ = 0.0;
  sSd_ = 0.0;
  seMean_ = 0.0;
}

RunStats::RunStats(const std::vector<double> & init)
  : RunStats() {
  update(init);
}


void RunStats::update(const double & add){
  if(n == 0){
    ++n;
    sMean_ = add;
    sSqMean_ = add*add;
    sVar_ = 0.0;
    sSd_ = 0.0;
    seMean_ = 0.0;
  }
  else{
    ++n;
    sMean_ = sMean_ + (add - sMean_)/double(n);
    sSqMean_ = sSqMean_ + (add*add - sSqMean_)/double(n);
    sVar_ = (sSqMean_ - sMean_*sMean_)*double(n)/double(n-1);
    sSd_ = std::sqrt(sVar_);
    seMean_ = sSd_/std::sqrt(double(n));
  }
}


void RunStats::update(const std::vector<double> & add){
  std::for_each(add.begin(),add.end(),
		[this](const double & x){
		  update(x);
		});
}


double RunStats::sMean() const {
  return sMean_;
}


double RunStats::sVar() const {
  return sVar_;
}


double RunStats::sSd() const {
  return sSd_;
}


double RunStats::seMean() const {
  return seMean_;
}
