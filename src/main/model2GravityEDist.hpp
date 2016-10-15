#ifndef MODEL_2_GRAVITY_E_DIST_HPP
#define MODEL_2_GRAVITY_E_DIST_HPP


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramBeta2.hpp"
#include "paramGravityEDist.hpp"
#include "paramTrt.hpp"
#include "mcmcGravity2.hpp"


class Model2GravityEDist : public ModelBase {
protected:
public:
    Model2GravityEDist();
    Model2GravityEDist(const FixedData & fD);
    Model2GravityEDist(const Model2GravityEDist & m);

    virtual Model2GravityEDist & operator=(const Model2GravityEDist & m);


    Gravity2Mcmc mcmc;
};


#endif
