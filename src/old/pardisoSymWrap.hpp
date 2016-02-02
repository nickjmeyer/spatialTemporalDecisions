#ifdef NJM_USE_MKL
#ifndef PARDISO_SYM_WRAP_HPP__
#define PARDISO_SYM_WRAP_HPP__

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>

#include "utilities.hpp"
#include "mkl_pardiso.h"
#include "mkl_types.h"

std::vector<double> pardisoSymWrap(const std::vector<int> & iaVec,
				   const std::vector<int> & jaVec,
				   const std::vector<double> & aVec,
				   const std::vector<double> & bVec);



void mat2Raw(const Eigen::SparseMatrix<double> & mat,
	     std::vector<int> & iaVec,
	     std::vector<int> & jaVec,
	     std::vector<double> & aVec);

Eigen::VectorXd pardisoSolve(const Eigen::SparseMatrix<double> & mat,
			     const Eigen::VectorXd & vec);


#endif
#endif
