#include "model2EdgeToEdge.hpp"

static std::vector<ParamBase *> genPars(){
    std::vector<ParamBase *> pars;
    pars.push_back(new ParamIntercept);
    pars.push_back(new ParamBeta2);
    pars.push_back(new ParamTrt);
    return pars;
}

Model2EdgeToEdge::Model2EdgeToEdge(const FixedData & fD)
    : ModelBase("2EdgeToEdge",genPars(),fD){
}


Model2EdgeToEdge::Model2EdgeToEdge()
    : ModelBase("2EdgeToEdge",genPars()){
}


Model2EdgeToEdge::Model2EdgeToEdge(const Model2EdgeToEdge & m){
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
    // numInfected = m.numInfected;
    // numNotInfec = m.numNotInfec;
    fitType = m.fitType;
    fixSample = m.fixSample;
    mcmc = m.mcmc;
    setEdgeToEdge(m.getEdgeToEdge());
}


Model2EdgeToEdge &
Model2EdgeToEdge::operator=(const Model2EdgeToEdge & m){
    if(this != & m){
        this->Model2EdgeToEdge::~Model2EdgeToEdge();
        new (this) Model2EdgeToEdge(m);
    }
    return *this;
}
