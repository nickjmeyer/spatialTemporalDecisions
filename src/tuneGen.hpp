#ifndef TUNE_GEN_HPP__
#define TUNE_GEN_HPP__


#include <iostream>
#include <numeric>
#include <eigen3/Eigen/Sparse>
#include <gsl/gsl_multimin.h>
#include "utilities.hpp"
#include "rand.hpp"
#include "data.hpp"
#include "system.hpp"
#include "modelParam.hpp"
#include "modelParamEbola.hpp"
#include "modelParamRange.hpp"
#include "modelParamCave.hpp"
#include "model.hpp"
#include "modelEbola.hpp"
#include "modelRange.hpp"
#include "modelCave.hpp"
#include "dataDepth.hpp"
#include "noTrtAgent.hpp"
#include "myopicAgent.hpp"
#include "rankAgentToy.hpp"
#include "runner.hpp"
#include "settings.hpp"
#include "m1SgdOptim.hpp"
#include "m1NmOptim.hpp"
#include "m2NmOptim.hpp"
#include "anchorMan.hpp"
#include "toyFeatures0.hpp"
#include "toyFeatures1.hpp"
#include "toyFeatures2.hpp"

typedef GravityModel GM;
typedef GravityParam GP;
typedef GM EM;
typedef GP EP;

typedef System<GM,GP,EM,EP> S;
typedef NoTrt<EM,EP> NT;
typedef MyopicAgent<EM,EP> MA;

typedef VanillaRunnerNS<S,NT> RN;
typedef VanillaRunnerNS<S,MA> RM;


double TuneGenNTObj(const gsl_vector * x, void * param);
double TuneGenMAObj(const gsl_vector * x, void * param);


#endif
