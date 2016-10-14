#ifndef TUNE_PARAM_HPP
#define TUNE_PARAM_HPP

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
