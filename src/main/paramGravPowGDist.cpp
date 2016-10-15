#include "paramGravPowGDist.hpp"



unsigned int ParamGravPowGDist::initParsSize(const FixedData & fD){
    return 3U;
}


std::vector<std::string> ParamGravPowGDist::initNames(){
    return {"alpha","power","gPow"};
}


std::vector<bool> ParamGravPowGDist::initToScale(){
    return {true,false,false};
}


void ParamGravPowGDist::initInternal(const FixedData & fD){
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


void ParamGravPowGDist::updateBefore(){
}


void ParamGravPowGDist::updateAfter(){
    double alpha = pars.at(0);
    double power = pars.at(1);
    double gPow = pars.at(2);
    int i,I = numNodes * numNodes;
    for(i = 0; i < I; ++i){
        grav.at(i) = alpha * std::pow(dist.at(i),gPow) / std::pow(cc.at(i),
                std::exp(power));
    }
}


void ParamGravPowGDist::setFill(std::vector<double> & probs,
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


void ParamGravPowGDist::setFill(std::vector<double> & probs,
        std::vector<double> & pcPartial,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    int i,I = numNodes*numNodes;
    std::vector<double>::iterator it0;
    std::vector<double>::const_iterator it1;
    const double alpha = pars.at(0);
    const double power = pars.at(1);
    const double gPow = pars.at(2);
    for(i = 0,
            it0 = probs.begin(),
            it1 = grav.begin(); i < I; ++it0,++it1,++i){ // not infected
        *it0 -= *it1;

        // the negative is because the term is subtracted in the model
        pcPartial.at(i*totNumPars + offset) =
            -std::pow(dist[i],gPow)/std::pow(cc[i],std::exp(power));
        // remember that this term is subtracted in the model
        // so the negatives cancel
        pcPartial.at(i*totNumPars + offset + 1) =
            alpha*std::pow(dist[i],gPow)*std::log(cc[i])*
            std::exp(power)/std::pow(cc[i],std::exp(power));

        pcPartial.at(i*totNumPars + offset + 2) =
            -alpha*std::pow(dist[i],gPow)*std::log(dist[i])/
                std::pow(cc[i],std::exp(power));
    }
}


void ParamGravPowGDist::modFill(std::vector<double> & probs,
				const SimData & sD,
				const TrtData & tD,
				const FixedData & fD,
				const DynamicData & dD){
}


void ParamGravPowGDist::modFill(std::vector<double> & probs,
        std::vector<double> & pcPartial,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
}


std::vector<double> ParamGravPowGDist::partial(const int notNode,
        const int infNode,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    double alpha = pars.at(0);
    double power = pars.at(1);
    double gPow = pars.at(2);
    std::vector<double> p;
    int ind = notNode*numNodes + infNode;
    // the negative is because the term is subtracted in the model
    p.push_back(-std::pow(dist[ind],gPow)/std::pow(cc[ind],std::exp(power)));
    // remember that this term is subtracted in the model
    // so the negatives cancel
    p.push_back(alpha*std::pow(dist[ind],gPow)*std::log(cc[ind])*
            std::exp(power)/std::pow(cc[ind],std::exp(power)));

    p.push_back(-alpha*std::pow(dist[ind],gPow)*std::log(dist[ind])/
            std::pow(cc[ind],std::exp(power)));
    return p;
}


std::vector<double> ParamGravPowGDist::partial2(const int notNode,
        const int infNode,
        const SimData & sD,
        const TrtData & tD,
        const FixedData & fD,
        const DynamicData & dD){
    double alpha = pars.at(0);
    double power = pars.at(1);
    double gPow = pars.at(2);
    std::vector<double> p;
    int ind = notNode*numNodes + infNode;

    p.push_back(0); // (0,0)
    // remember that this term is subtracted in the model
    // so the negatives cancel
    double val0 = std::pow(dist[ind],gPow)*std::log(cc[ind])*
        std::exp(power)/std::pow(cc[ind],std::exp(power));
    p.push_back(val0); // (0,1)

    double val1 = -std::log(dist[ind])*std::pow(dist[ind],gPow)/
        std::pow(cc[ind],std::exp(power));
    p.push_back(val1); // (0,2)

    p.push_back(val0);// (1,0)

    //(1,1)
    p.push_back(- alpha * val0 * (std::exp(power)*std::log(cc[ind]) - 1.0));

    double val3 = alpha * val0 * std::log(dist[ind]);
    p.push_back(val3); // (1,2)

    p.push_back(val1); // (2,0)

    p.push_back(val3); // (2,1)

    double val4 = alpha * val1 * std::log(dist[ind]);
    p.push_back(val4); // (2,2)



    return p;
}
