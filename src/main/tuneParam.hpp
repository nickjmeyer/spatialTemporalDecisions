#ifndef TUNE_PARAM_HPP__
#define TUNE_PARAM_HPP__

#include <vector>


class TuneParam {
public:
  virtual std::vector<double> getPar() const = 0;
  virtual void putPar(const std::vector<double> & par) = 0;


  bool getEdgeToEdge() const;

  void setEdgeToEdge(const bool edgeToEdge);

private:
  bool edgeToEdge;
};


#endif
