#ifndef MODEL_2_EDGE_TO_EDGE_HPP
#define MODEL_2_EDGE_TO_EDGE_HPP


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
#include "mcmcEdgeToEdge2.hpp"


class Model2EdgeToEdge : public ModelBase {
protected:
public:
    Model2EdgeToEdge();
    Model2EdgeToEdge(const FixedData & fD);
    Model2EdgeToEdge(const Model2EdgeToEdge & m);

    virtual Model2EdgeToEdge & operator=(const Model2EdgeToEdge & m);

    EdgeToEdge2Mcmc mcmc;
};


#endif
