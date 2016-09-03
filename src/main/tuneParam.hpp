#ifndef TUNE_PARAM_HPP__
#define TUNE_PARAM_HPP__

#include <vector>


class TuneParam {
  virtual std::vector<double> getPar() const = 0;
  virtual void putPar(const std::vector<double> & par) = 0;
};


#endif
