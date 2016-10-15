#include "model2GravityEDist.hpp"

static std::vector<ParamBase *> genPars(){
    std::vector<ParamBase *> pars;
    pars.push_back(new ParamIntercept);
    pars.push_back(new ParamBeta2);
    pars.push_back(new ParamGravityEDist);
    pars.push_back(new ParamTrt);
    return pars;
}

Model2GravityEDist::Model2GravityEDist(const FixedData & fD)
    : ModelBase("2GravityEDist",genPars(),fD){
}


Model2GravityEDist::Model2GravityEDist()
    : ModelBase("2GravityEDist",genPars()){
}


Model2GravityEDist::Model2GravityEDist(const Model2GravityEDist & m){
    int i, parsSize = m.pars.size();
    pars.clear();
    for(i = 0; i < parsSize; ++i)
        pars.push_back(m.pars.at(i)->clone());

    name = m.name;
    numPars = m.numPars;
    set = m.set;
    probs = m.probs;
    expitInfProbs = m.expitInfProbs;
    expitRevProbs = m.expitRevProbs;
    quick = m.quick;
    pcPartial = m.pcPartial;
    meanHit = m.meanHit;
    varHit = m.varHit;
    ready = m.ready;
    numInfected = m.numInfected;
    numNotInfec = m.numNotInfec;
    fitType = m.fitType;
    fixSample = m.fixSample;
    mcmc = m.mcmc;
    setEdgeToEdge(m.getEdgeToEdge());
}


Model2GravityEDist &
Model2GravityEDist::operator=(const Model2GravityEDist & m){
    if(this != & m){
        this->Model2GravityEDist::~Model2GravityEDist();
        new (this) Model2GravityEDist(m);
    }
    return *this;
}
