#ifndef SYSTEM_HPP
#define SYSTEM_HPP


#include <vector>
#include <string>
#include "utilities.hpp"
#include "timer.hpp"
#include "data.hpp"
#include "model.hpp"
#include "modelGravityGDist.hpp"
#include "model2GravityGDist.hpp"
#include "model2GravityEDist.hpp"
#include "model2GPowGDist.hpp"
#include "model2EdgeToEdge.hpp"
#include "modelGDist.hpp"
#include "modelEDist.hpp"
#include "modelIntercept.hpp"
#include "rand.hpp"
#include "calcCentrality.hpp"

// used in precomputing data
struct ExpDistData {
    std::vector<double> dist;
    int cutoff;
    double proportion;
};

double expDistEval(double c, void * params);
double expDistGrad(double c, void * params);
void expDistEvalGrad(double c, void *params, double *eval, double *grad);

template <class MG, class ME>
class System {
public:

    System();
    System(const SimData & sD, const TrtData & tD,
            const FixedData & fD, const DynamicData & dD,
            const MG & modelGen, const ME & modelEst);
    System(const std::string file);

    int specialInit;

    SimData sD;
    TrtData tD;
    FixedData fD;
    DynamicData dD;

    SimData sD_r;
    TrtData tD_r;
    DynamicData dD_r;

    MG modelGen;
    ME modelEst;

    MG modelGen_r;
    ME modelEst_r;


    virtual void reset(const std::vector<int> & ind);

    virtual void revert();

    virtual void checkPoint();

    virtual void initialize();

    virtual void setEdgeToEdge(const bool edgeToEdge);

    virtual void preCompData();

    virtual void nextPoint();
    virtual void nextPoint(const std::vector<double> & infProbs);

    virtual void updateStatus();

    virtual double value();
};





#endif
