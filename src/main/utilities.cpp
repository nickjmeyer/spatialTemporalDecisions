#include "utilities.hpp"

double njm::expit(const double n){
  return 1.0 - 1.0/(1.0 + std::exp(n));
}


double njm::l2norm(const std::vector<double> & v0,
	      const std::vector<double> & v1){
  double norm=0;
  int i,size = std::min(v0.size(),v1.size());
  for(i=0; i<size; i++)
    norm += (v0.at(i) - v1.at(i))*(v0.at(i) - v1.at(i));
  return std::sqrt(norm);
}

double njm::l2norm(std::vector<double> & v){
  double norm=0;
  int i,size = v.size();
  for(i=0; i<size; i++)
    norm += v.at(i)*v.at(i);
  norm = std::sqrt(norm);
  if(norm > 0)
    for(i=0; i<size; i++)
      v.at(i)/=norm;
  return norm;
}

double njm::sampVar(const std::vector<double> & v){
  std::vector<double>::const_iterator it,end;
  end = v.end();
  double mean = 0, meanSq = 0;
  for (it = v.begin(); it != end; ++it) {
    const double val = *it;
    mean += val;
    meanSq += val*val;
  }
  const int cnt = v.size();
  mean /= cnt;
  meanSq /= cnt;
  return (double(cnt) / double(cnt - 1)) * (meanSq - mean*mean);
}
