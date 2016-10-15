#include "bayesPsamplesBR.hpp"


std::vector<double> getStats(const std::vector<std::vector<int> > & h,
        const SimData & sD,
        const FixedData & fD){
    std::vector<double> stats;
    // total infected
    stats.push_back(sD.numInfected);

    int i,j,sum = 0;
    for(i = 1; i < (int)h.size(); ++i){
        sum = 0;
        for(j = 0; j < fD.numNodes; ++j)
            if(h.at(i).at(j) >= 2 && h.at(i-1).at(j) < 2)
                ++sum;

        stats.push_back(sum);
    }


    // average year of infection
    sum = 0;
    for(i = 0; i < (int)h.size(); ++i){
        for(j = 0; j < fD.numNodes; ++j){
            if((i == 0 && h.at(i).at(j) >= 2) ||
                    (i > 0 && h.at(i).at(j) >= 2 && h.at(i-1).at(j) < 2))
                sum += i;
        }
    }
    stats.push_back(((double)sum)/((double)sD.numInfected));


    // average lat long for infected
    double cLong,cLat;
    cLong = cLat = 0;
    for(i = 0 ; i < sD.numInfected; ++i){
        cLong += fD.centroidsLong.at(sD.infected.at(i));
        cLat += fD.centroidsLat.at(sD.infected.at(i));
    }
    stats.push_back(cLong/((double)sD.numInfected));
    stats.push_back(cLat/((double)sD.numInfected));

    // average spread distance from start for infected
    double sDist;
    int start = -1,numStart = 0;
    for(i = 0; i < fD.numNodes; ++i){
        if(h.at(0).at(i) >= 2){
            if(numStart == 0)
                start = i;
            ++numStart;
        }
    }
    if(numStart > 1)
        std::cout << "warning...multiple starting locations..."
                  << std::endl;
    sDist = 0;
    for(i = 0; i < sD.numInfected; ++i){
        sDist += fD.gDist.at(sD.infected.at(i)*fD.numNodes + start);
    }
    sDist/=(double)sD.numInfected;
    stats.push_back(sDist);

    // min lat long for infected
    cLong = std::numeric_limits<double>::max();
    cLat = std::numeric_limits<double>::max();
    for(i = 0; i < sD.numInfected; ++i){
        if(fD.centroidsLong.at(sD.infected.at(i)) < cLong)
            cLong = fD.centroidsLong.at(sD.infected.at(i));
        if(fD.centroidsLat.at(sD.infected.at(i)) < cLat)
            cLat = fD.centroidsLat.at(sD.infected.at(i));
    }
    stats.push_back(cLong);
    stats.push_back(cLat);


    // max lat long for infected
    cLong = std::numeric_limits<double>::lowest();
    cLat = std::numeric_limits<double>::lowest();
    for(i = 0; i < sD.numInfected; ++i){
        if(fD.centroidsLong.at(sD.infected.at(i)) > cLong)
            cLong = fD.centroidsLong.at(sD.infected.at(i));
        if(fD.centroidsLat.at(sD.infected.at(i)) > cLat)
            cLat = fD.centroidsLat.at(sD.infected.at(i));
    }
    stats.push_back(cLong);
    stats.push_back(cLat);


    // max dist from starting location
    sDist = std::numeric_limits<double>::lowest();
    for(i = 0; i < sD.numInfected; ++i)
        if(fD.gDist.at(sD.infected.at(i)*fD.numNodes + start) > sDist)
            sDist = fD.gDist.at(sD.infected.at(i)*fD.numNodes + start);
    stats.push_back(sDist);


    return stats;
}




int main(int argc, char ** argv){
    njm::sett.set(argc,argv);


    typedef Model2GravityGDist MG;

    typedef MG ME;

    typedef System<MG,ME> S;


    std::string sampsFile = "Param2GravityGDistXOXOXO/samples.txt";
    std::vector<double> samples;
    njm::fromFile(samples,njm::sett.srcExt(sampsFile));


    std::vector<std::string> names = {"n_inf","n_inf_2007","n_inf_2008",
                                      "n_inf_2009","n_inf_2010","n_inf_2011",
                                      "n_inf_2012","n_inf_2013","mean_year",
                                      "mean_long","mean_lat",
                                      "mean_dist_from_start",
                                      "min_long","min_lat",
                                      "max_long","max_lat",
                                      "max_dist_from_start"};

    std::vector<std::vector<double> > stats;

    njm::toFile(names,njm::sett.datExt("sampStats_",".txt"),
            std::ios_base::out);

    int numStats = 10000;


    S s;
    Starts starts("startingLocations.txt");
    int r,t,R,T;
    std::vector<std::vector<int> > h;

    int numSamps = samples.size()/s.modelGen_r.numPars;
    assert(numSamps*((int)s.modelGen_r.numPars) == (int)samples.size());

    R = numStats;
    T = s.fD.trtStart;
    for(r = 0; r < R; ++r){
        printf("% 8d\r",r);
        fflush(stdout);

        std::vector<double> par;
        int samp = njm::runifInterv(0,numSamps);
        int i;
        for(i = 0; i < (int)s.modelGen_r.numPars; ++i)
            par.push_back(samples.at(samp*s.modelGen_r.numPars + i));

        s.modelGen_r.putPar(par.begin());
        s.modelEst_r.putPar(par.begin());

        s.reset(starts[r]);

        for(t = 0; t < T; ++t)
            s.nextPoint();
        h = s.sD.history;
        h.push_back(s.sD.status);

        stats.push_back(getStats(h,s.sD,s.fD));
    }

    printf("done               \n");

    njm::toFile(njm::toString(stats,"\n",""),
            njm::sett.datExt("sampStats_",".txt"));


    return 0;
}
