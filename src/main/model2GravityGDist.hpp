#ifndef MODEL_2_GRAVITY_G_DIST_HPP
#define MODEL_2_GRAVITY_G_DIST_HPP


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramBeta2.hpp"
#include "paramGravityGDist.hpp"
#include "paramTrt.hpp"

class Model2GravityGDist : public ModelBase {
protected:
public:
    Model2GravityGDist();
    Model2GravityGDist(const FixedData & fD);
    Model2GravityGDist(const Model2GravityGDist & m);

    virtual Model2GravityGDist & operator=(const Model2GravityGDist & m);
};


#endif
