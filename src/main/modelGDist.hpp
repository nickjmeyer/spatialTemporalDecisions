#ifndef MODEL_G_DIST_HPP
#define MODEL_G_DIST_HPP


#include <armadillo>
#include <cmath>
#include <gsl/gsl_multimin.h>
#include "data.hpp"
#include "settings.hpp"
#include "model.hpp"
#include "paramIntercept.hpp"
#include "paramGDist.hpp"
#include "paramTrt.hpp"


class ModelGDist : public ModelBase {
protected:
public:
    ModelGDist();
    ModelGDist(const FixedData & fD);
    ModelGDist(const ModelGDist & m);

    virtual ModelGDist & operator=(const ModelGDist & m);
};


#endif
