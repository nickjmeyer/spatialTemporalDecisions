#include "paramBeta.hpp"


unsigned int ParamBeta::initParsSize(const FixedData & fD){
    return fD.numCovar; // unsigned int literal
}


std::vector<std::string> ParamBeta::initNames(){
    std::vector<std::string> str;
    unsigned int i;
    for(i = 0; i < parsSize; ++i)
        str.push_back("beta"+njm::toString(i,"",0,0));
    return str;
}


void ParamBeta::initInternal(const FixedData & fD){
    numNodes = fD.numNodes;
    covar = fD.covar;
    covarBeta = std::vector<double>(numNodes,0);
}


void ParamBeta::updateBefore(){
}


void ParamBeta::updateAfter(){
    int i,j,k,J = parsSize;
    for(i = 0,k = 0; i < numNodes; ++i){
        covarBeta.at(i) = 0;
        for(j = 0; j < J; ++j, ++k){
            covarBeta.at(i) += covar.at(k) * pars.at(j);
        }
    }
}


void ParamBeta::setFill(std::vector<double> & probs,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    int i,j;
    std::vector<double>::iterator it;
    for(i = 0, it = probs.begin(); i < numNodes; ++i){ // not infected
        for(j = 0; j < numNodes; ++j, ++it){ // infected
            *it += covarBeta.at(i);
        }
    }
}


void ParamBeta::setFill(std::vector<double> & probs,
        std::vector<double> & pcPartial,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    int i,j,k;
    std::vector<double>::iterator it;
    for(i = 0, it = probs.begin(); i < numNodes; ++i){ // not infected
        for(j = 0; j < numNodes; ++j, ++it){ // infected
            *it += covarBeta.at(i);

            // set partial derivative
            for(k = 0; k < int(parsSize); ++k) {
                pcPartial.at(i*numNodes*totNumPars + j*totNumPars
                        + offset + k) = covar.at(i*parsSize + k);
            }
        }
    }
}


void ParamBeta::modFill(std::vector<double> & probs,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
}


void ParamBeta::modFill(std::vector<double> & probs,
        std::vector<double> & pcPartial,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
}


std::vector<double> ParamBeta::partial(const int notNode,
        const int infNode,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    std::vector<double> p;
    int i,ind = notNode*parsSize;
    for(i = 0; i < int(parsSize); ++i){
        p.push_back(covar.at(ind + i));
    }
    return p;
}
