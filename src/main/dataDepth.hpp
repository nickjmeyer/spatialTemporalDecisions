#ifndef DATA_DEPTH_HPP
#define DATA_DEPTH_HPP


#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include "sortMerge.hpp"

double halfPlaneDepth(const double U, const double V, const int N,
        const std::vector<double> & X,
        const std::vector<double> & Y);


#endif
