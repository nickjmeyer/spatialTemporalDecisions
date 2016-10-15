#include "modelGravityGDist.hpp"

static std::vector<ParamBase *> genPars(){
    std::vector<ParamBase *> pars;
    pars.push_back(new ParamIntercept);
    pars.push_back(new ParamBeta);
    pars.push_back(new ParamGravityGDist);
    pars.push_back(new ParamTrt);
    return pars;
}

ModelGravityGDist::ModelGravityGDist(const FixedData & fD)
    : ModelBase("GravityGDist",genPars(),fD){
}


ModelGravityGDist::ModelGravityGDist()
    : ModelBase("GravityGDist",genPars()){
}


ModelGravityGDist::ModelGravityGDist(const ModelGravityGDist & m){
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
    mcmc = m.mcmc;
    fixSample = m.fixSample;
    setEdgeToEdge(m.getEdgeToEdge());
}


ModelGravityGDist & ModelGravityGDist::operator=(const ModelGravityGDist & m){
    if(this != & m){
        this->ModelGravityGDist::~ModelGravityGDist();
        new (this) ModelGravityGDist(m);
    }
    return *this;
}



double ModelGravityGDist::tuneTrt(const FixedData & fD){
    int i,j;
    double avgCaves = 0.0;
    for(i = 0; i < fD.numNodes; i++)
        avgCaves += fD.caves.at(i);
    avgCaves /= double(fD.numNodes);

    double minDist = std::numeric_limits<double>::max();
    for(i = 0; i < fD.numNodes; i++)
        for(j = (i+1); j < fD.numNodes; j++)
            if(minDist > fD.gDist.at(i*fD.numNodes + j))
                minDist = fD.gDist.at(i*fD.numNodes + j);

    std::vector<double> vals;
    vals = pars.at(0)->getPar();
    double base = vals.at(0);
    vals = pars.at(2)->getPar();
    double alpha = vals.at(0);
    double power = vals.at(1);
    base -= alpha * minDist/std::pow(avgCaves*avgCaves,std::exp(power));

    return -(std::log(0.005) - base)/2.0;
}
