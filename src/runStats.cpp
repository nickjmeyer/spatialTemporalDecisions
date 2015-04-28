#include "runStats.hpp"


RunStats::RunStats(const std::vector<double> & init){
  n = 0U;
  vals.clear();
  smean_ = 0.0;
  svar_ = 0.0;
  ssd_ = 0.0;
}


void RunStats::update(const double & add){
  if(n == 0){
    ++n;
    smean_ = add;
    svar_ = 0.0;
    ssd_ = 0.0;
    seMean_ = 0.0;
  }
  else{
    ++n;
    svar_ = svar_*double(n-2)/double(n-1)
      + (add - smean_)*(add - smean_)/double(n);
    smean_ = smean_ + (add - smean_)/double(n);
    ssd_ = std::sqrt(svar_);
    seMean_ = ssd_/std::sqrt(double(n));
  }
}


void RunStats::update(const std::vector<double> & add){
  std::for_each(add.begin(),add.end(),
		[this](const double & x){
		  update(x);
		});
}


double RunStats::smean() const {
  return smean_;
}


double RunStats::svar() const {
  return svar_;
}


double RunStats::ssd() const {
  return ssd_;
}


double RunStats::seMean() const {
  return seMean_;
}
