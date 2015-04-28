#include "runStats.hpp"


RunStats::RunStats(const std::vector<double> & init){
  n = 0U;
  vals.clear();
  smean = 0.0;
  svar = 0.0;
  ssd = 0.0;
}


void RunStats::update(const double & add){
  if(n == 0){
    smean = add;
    svar = 0.0;
    ssd = 0.0;
  }
  else{
    ++n;
    svar = (n-2)svar/(n-1) + (add - smean)*(add-smean)/n;
  }
}

