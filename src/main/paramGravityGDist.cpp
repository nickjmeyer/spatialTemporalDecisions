#include "paramGravityGDist.hpp"



unsigned int ParamGravityGDist::initParsSize(const FixedData & fD){
    return 2U;
}


std::vector<std::string> ParamGravityGDist::initNames(){
    return {"alpha","power"};
}


std::vector<bool> ParamGravityGDist::initToScale(){
    return {true,false};
}


void ParamGravityGDist::initInternal(const FixedData & fD){
    numNodes = fD.numNodes;
    grav = std::vector<double>(numNodes*numNodes,0.0);

    dist = fD.gDist;

    cc.clear();
    cc.reserve(numNodes*numNodes);
    int i,j;
    for(i = 0; i < numNodes; ++i){
        for(j = 0; j < numNodes; ++j){
            cc.push_back(fD.caves.at(i)*fD.caves.at(j));
        }
    }
}


void ParamGravityGDist::updateBefore(){
}


void ParamGravityGDist::updateAfter(){
    double alpha = pars.at(0);
    double power = pars.at(1);
    int i,I = numNodes * numNodes;
    for(i = 0; i < I; ++i){
        grav.at(i) = alpha * dist.at(i) / std::pow(cc.at(i),std::exp(power));
    }
}


void ParamGravityGDist::setFill(std::vector<double> & probs,
				const SimData & sD,
				const TrtData & tD,
				const FixedData & fD,
				const DynamicData & dD){
    int i,I = numNodes*numNodes;
    std::vector<double>::iterator it0;
    std::vector<double>::const_iterator it1;
    for(i = 0,
            it0 = probs.begin(),
            it1 = grav.begin(); i < I; ++it0,++it1,++i){ // not infected
        *it0 -= *it1;
    }
}


void ParamGravityGDist::modFill(std::vector<double> & probs,
				const SimData & sD,
				const TrtData & tD,
				const FixedData & fD,
				const DynamicData & dD){
}


std::vector<double> ParamGravityGDist::partial(const int notNode,
        const int infNode,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    double alpha = pars.at(0);
    double power = pars.at(1);
    std::vector<double> p;
    int ind = notNode*numNodes + infNode;
    // the negative is because the term is subtracted in the model
    p.push_back(-dist[ind]/std::pow(cc[ind],std::exp(power)));
    // remember that this term is subtracted in the model
    // so the negatives cancel
    p.push_back(alpha*dist[ind]*std::log(cc[ind])*std::exp(power)/
            std::pow(cc[ind],std::exp(power)));
    return p;
}



std::vector<double> ParamGravityGDist::partial2(const int notNode,
        const int infNode,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    double alpha = pars.at(0);
    double power = pars.at(1);
    std::vector<double> p;
    int ind = notNode*numNodes + infNode;

    p.push_back(0);
    // remember that this term is subtracted in the model
    // so the negatives cancel
    double val = dist[ind]*std::log(cc[ind])*std::exp(power)/
        std::pow(cc[ind],std::exp(power));

    p.push_back(val);
    p.push_back(val);

    p.push_back(- alpha * val * (std::exp(power)*std::log(cc[ind]) - 1.0));
    return p;
}
