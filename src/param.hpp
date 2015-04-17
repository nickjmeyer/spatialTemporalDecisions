#ifndef PARAM_HPP__
#define PARAM_HPP__

#include <vector>
#include <algorithm>
#include <eigen3/Eigen/Eigen>
#include "data.hpp"

class ParamBase {
 protected:
  
  std::vector<double> pars;
  std::vector<double>::iterator beg,end;
  unsigned int parsSize;


  virtual unsigned int initParsSize(const FixedData & fD) = 0;
  
  // initialize the internal book-keeping data
  virtual void initInternal(const FixedData & fD) = 0;
  
  // update the internal book-keeping data
  // called inside of putPar() to guarantee up-to-date book-keeping data
  virtual void updateBefore() = 0;
  virtual void updateAfter() = 0;

  // initializes pars = {0,...}, beg = pars.begin(), end = pars.end()
  // calls initInternal, initParsSize
  // should be called from derived constructor
  virtual void init(const FixedData & fD);
  
 public:
  virtual ~ParamBase() { };
  
  // retreive pars
  std::vector<double> getPar() const;
  
  // change pars
  // requires that newParInt has at least parsSize elements left
  // requires that beg,end point to begin() and end() of pars
  std::vector<double>::iterator putPar(std::vector<double>::iterator newParIt);

  // set the probs matrix
  virtual void setFill(std::vector<double> & probs,
		       const SimData & sD,
		       const TrtData & tD,
		       const FixedData & fD,
		       const DynamicData & dD) const = 0;

  // update the probs matrix
  virtual void updateFill(std::vector<double> & probs,
			  const SimData & sD,
			  const TrtData & tD,
			  const FixedData & fD,
			  const DynamicData & dD) const = 0;
};


#endif
