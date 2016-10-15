#include "obsDataStats.hpp"


int main(int argc, char ** argv){
    njm::sett.set(argc,argv);

    typedef ModelTimeExpCaves MG;
    typedef MG ME;
    typedef System<MG,ME> S;
    typedef WnsFeatures1<ME> F;
    typedef RankAgent<F,ME> RA;
    typedef M1SpOptim<S,RA,ME> SPO;

    S s("obsData.txt");
    // S s;
    // Starts starts("startingLocations.txt");
    // s.reset(starts[0]);
    s.modelEst_r = s.modelGen_r;
    s.revert();

    F f;
    RA ra;
    SPO spo;
    spo.tp.tune = 0;

    // optimize to get weights
    std::cout << "optimizing..." << std::flush;
    spo.optim(s,ra);
    std::cout << "done" << std::endl;

    // apply trt
    std::cout << "trt..." << std::flush;
    ra.applyTrt(s.sD,s.tD,s.fD,s.dD,s.modelEst);
    std::cout << "done" << std::endl;

    // get features
    std::cout << "features..." << std::flush;
    f.preCompData(s.sD,s.tD,s.fD,s.dD,
            s.modelEst);
    f.getFeatures(s.sD,s.tD,s.fD,s.dD,
            s.modelEst);
    std::cout << "done" << std::endl;

    // calculate priority scores
    std::cout << "scores..." << std::flush;
    arma::colvec notPs,infPs;
    notPs = f.notFeat*ra.tp.weights;
    infPs = f.infFeat*ra.tp.weights;

    std::vector<double> ps(s.fD.numNodes,0);
    int i;
    for(i = 0; i < s.sD.numNotInfec; ++i){
        ps.at(s.sD.notInfec.at(i)) = notPs(i);
    }
    for(i = 0; i < s.sD.numInfected; ++i){
        ps.at(s.sD.infected.at(i)) = infPs(i);
    }
    std::cout << "done" << std::endl;

    // save to file
    njm::toFile(njm::toString(ps,"\n",""),njm::sett.srcExt("obsDataPs.txt"),
            std::ios_base::out);



    njm::sett.clean();
    return 0;
}
