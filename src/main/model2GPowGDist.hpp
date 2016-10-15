#ifndef MODEL_2_G_POW_G_DIST_HPP
#define MODEL_2_G_POW_G_DIST_HPP


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramBeta2.hpp"
#include "paramGravPowGDist.hpp"
#include "paramTrt.hpp"


class Model2GPowGDist : public ModelBase {
protected:
public:
    Model2GPowGDist();
    Model2GPowGDist(const FixedData & fD);
    Model2GPowGDist(const Model2GPowGDist & m);

    virtual Model2GPowGDist & operator=(const Model2GPowGDist & m);
};


#endif
