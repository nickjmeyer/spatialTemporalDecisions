#ifndef BAYES_P_SAMPLES_BR_HPP__
#define BAYES_P_SAMPLES_BR_HPP__


#include <iostream>
#include <vector>
#include <string>
#include <assert.h>
#include "utilities.hpp"
#include "rand.hpp"
#include "data.hpp"
#include "system.hpp"
#include "starts.hpp"



std::vector<double> getStats(const std::vector<std::vector<int> > & h,
			     const SimData & sD,
			     const FixedData & fD);



#endif
