#include "paramBeta2.hpp"
#include "timer.hpp"


unsigned int ParamBeta2::initParsSize(const FixedData & fD){
    return fD.numCovar*2; // unsigned int literal
}


std::vector<std::string> ParamBeta2::initNames(){
    std::vector<std::string> str;
    unsigned int i;
    for(i = 0; i < parsSize; ++i)
        str.push_back("beta"+njm::toString(i,"",0,0));
    return str;
}


void ParamBeta2::initInternal(const FixedData & fD){
    numNodes = fD.numNodes;
    covar = fD.covar;
    numCovar = fD.numCovar;
    covarBeta = std::vector<double>(numNodes,0);
    covarBetaInf = std::vector<double>(numNodes,0);
}


void ParamBeta2::updateBefore(){
}


void ParamBeta2::updateAfter(){
    int i,j,k,J = numCovar;
    for(i = 0,k = 0; i < numNodes; ++i){
        covarBeta.at(i) = 0;
        covarBetaInf.at(i) = 0;
        for(j = 0; j < J; ++j, ++k){
            covarBeta.at(i) += covar.at(k) * pars.at(j);
            covarBetaInf.at(i) += covar.at(k) * pars.at(j+numCovar);
        }
    }
}


void ParamBeta2::setFill(std::vector<double> & probs,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    int i,j;
    std::vector<double>::iterator it;
    for(i = 0, it = probs.begin(); i < numNodes; ++i){ // not infected
        for(j = 0; j < numNodes; ++j, ++it){ // infected
            *it += covarBeta.at(i) + covarBetaInf.at(j);
        }
    }
}


void ParamBeta2::setFill(std::vector<double> & probs,
        std::vector<double> & pcPartial,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    int i,j,k;
    std::vector<double>::iterator it;
    for(i = 0, it = probs.begin(); i < numNodes; ++i){ // not infected
        for(j = 0; j < numNodes; ++j, ++it){ // infected
            *it += covarBeta.at(i) + covarBetaInf.at(j);

            // set partial derivative
            for(k = 0; k < numCovar; ++k) {
                pcPartial.at(i*numNodes*totNumPars + j*totNumPars
                        + offset + k) = covar.at(i*numCovar + k);
            }
            for(k = 0; k < numCovar; ++k) {
                pcPartial.at(i*numNodes*totNumPars + j*totNumPars
                        + offset + numCovar + k) = covar.at(j*numCovar + k);
            }
        }
    }
}


void ParamBeta2::modFill(std::vector<double> & probs,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
}


void ParamBeta2::modFill(std::vector<double> & probs,
        std::vector<double> & pcPartial,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
}


std::vector<double> ParamBeta2::partial(const int notNode,
        const int infNode,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    // njm::timer.start("partial_beta2");
    std::vector<double> p;

    int i,ind = notNode*numCovar,indInf = infNode*numCovar;
    for(i = 0; i < numCovar; ++i){
        p.push_back(covar.at(ind + i));
    }
    for(i = 0; i < numCovar; ++i){
        p.push_back(covar.at(indInf + i));
    }
    // njm::timer.stop("partial_beta2");
    return p;
}
