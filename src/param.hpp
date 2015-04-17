#ifndef PARAM_HPP__
#define PARAM_HPP__

#include <vector>
#include <algorithm>
#include <eigen3/Eigen/Eigen>
#include "data.hpp"

class BaseParam {
 protected:
  
  std::vector<double> pars;
  std::vector<double>::iterator beg,end;
  unsigned int parsSize;


  virtual unsigned int initParsSize(const FixedData & fD) = 0;
  
  // initialize the internal book-keeping data
  virtual void initInternal(const FixeData & fD) = 0;
  
  // update the internal book-keeping data
  // called inside of putPar() to guarantee up-to-date book-keeping data
  virtual void updateInternal() = 0;
  
 public:
  // initializes pars = {0,...}, beg = pars.begin(), end = pars.end()
  // calls initInternal, initParsSize
  BaseParam(const FixedData & fD);
  
  // retreive pars
  std::vector<double> getPar() const;
  
  // change pars
  // requires that newParInt has at least parsSize elements left
  // requires that beg,end point to begin() and end() of pars
  void putPar(std::vector<double>::iterator & newParIt);

  // modify the probs matrix
  virtual void updateFill(Eigen::MatrixXd & probs) = 0;
}


#endif
