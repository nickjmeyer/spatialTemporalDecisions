#include "bayesP.hpp"
#include <gflags/gflags.h>
#include <glog/logging.h>

using namespace google;
using namespace gflags;


DEFINE_string(srcDir,"","Path to source directory");
DEFINE_bool(dryRun,false,"Do not execute main");
DEFINE_bool(save,false,"Save parameter estimates for simluation");

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
        sDist += fD.eDist.at(sD.infected.at(i)*fD.numNodes + start);
    }
    sDist/=(double)sD.numInfected;
    stats.push_back(sDist*fD.eDistSd);

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
        if(fD.eDist.at(sD.infected.at(i)*fD.numNodes + start) > sDist)
            sDist = fD.eDist.at(sD.infected.at(i)*fD.numNodes + start);
    stats.push_back(sDist*fD.eDistSd);


    return stats;
}


std::vector<double> getStatsOos(
        const std::vector<std::vector<int> > & baseH,
        const std::vector<std::vector<int> > & newH,
        const FixedData & fD){
    std::vector<double> stats;
    // total infected
    int i,cnt = 0;
    for (i = 0; i < fD.numNodes; ++i) {
        if (newH.at(newH.size()-1).at(i) >= 2)
            ++cnt;
    }
    stats.push_back(cnt);

    int j,sum = 0;
    for(i = 0; i < (int)newH.size(); ++i){
        sum = 0;
        for(j = 0; j < fD.numNodes; ++j) {
            if(i == 0) {
                if(newH.at(i).at(j) >= 2 && baseH.at(baseH.size()-1).at(j) < 2)
                    ++sum;
            }
            else {
                if(newH.at(i).at(j) >= 2 && newH.at(i-1).at(j) < 2)
                    ++sum;
            }
        }

        stats.push_back(sum);
    }


    // average year of new infections
    sum = 0;
    cnt = 0;
    for(i = 0; i < (int)newH.size(); ++i){
        for(j = 0; j < fD.numNodes; ++j){
            if(i == 0) {
                if(newH.at(i).at(j) >= 2 && baseH.at(baseH.size()-1).at(j) < 2){
                    sum += baseH.size() + i;
                    ++cnt;
                }
            }
            else {
                if(newH.at(i).at(j) >= 2 && newH.at(i-1).at(j) < 2){
                    sum += baseH.size() + i;
                    ++cnt;
                }
            }
        }
    }
    stats.push_back(((double)sum)/((double)cnt));


    // average lat long for infected
    double cLong,cLat;
    cLong = cLat = 0;
    cnt = 0;
    for(i = 0 ; i < fD.numNodes; ++i){
        if(newH.at(newH.size()-1).at(i) >= 2 &&
                baseH.at(baseH.size()-1).at(i) < 2) {
            cLong += fD.centroidsLong.at(i);
            cLat += fD.centroidsLat.at(i);
            ++cnt;
        }
    }
    stats.push_back(cLong/((double)cnt));
    stats.push_back(cLat/((double)cnt));

    // average spread distance from start for infected
    std::vector<double> distFromStart;
    for (i = 0; i < fD.numNodes; ++i) {
        if (newH.at(newH.size()-1).at(i) >= 2 &&
                baseH.at(baseH.size()-1).at(i) < 2) {
            double minDistToStat = std::numeric_limits<double>::max();
            for (j = 0; j < fD.numNodes; ++j) {
                if(newH.at(newH.size()-1).at(j) >= 2 &&
                        baseH.at(baseH.size()-1).at(j) >= 2) {
                    minDistToStat = std::min(minDistToStat,
                            fD.eDist.at(i*fD.numNodes + j));
                }
            }
            distFromStart.push_back(minDistToStat);
        }
    }
    CHECK_EQ(distFromStart.size(),stats.at(1) + stats.at(2));

    double sumDist = 0;
    sumDist = std::accumulate(distFromStart.begin(),distFromStart.end(),0.);
    sumDist/=(double)distFromStart.size();
    stats.push_back(sumDist*fD.eDistSd);


    // min lat long for infected
    cLong = std::numeric_limits<double>::max();
    cLat = std::numeric_limits<double>::max();
    for(i = 0; i < fD.numNodes; ++i){
        if(newH.at(newH.size()-1).at(i) >= 2 &&
                baseH.at(baseH.size()-1).at(i) < 2) {
            if(fD.centroidsLong.at(i) < cLong)
                cLong = fD.centroidsLong.at(i);
            if(fD.centroidsLat.at(i) < cLat)
                cLat = fD.centroidsLat.at(i);
        }
    }
    stats.push_back(cLong);
    stats.push_back(cLat);


    // max lat long for infected
    cLong = std::numeric_limits<double>::lowest();
    cLat = std::numeric_limits<double>::lowest();
    for(i = 0; i < fD.numNodes; ++i){
        if(newH.at(newH.size()-1).at(i) >= 2 &&
                baseH.at(baseH.size()-1).at(i) < 2) {
            if(fD.centroidsLong.at(i) > cLong)
                cLong = fD.centroidsLong.at(i);
            if(fD.centroidsLat.at(i) > cLat)
                cLat = fD.centroidsLat.at(i);
        }
    }
    stats.push_back(cLong);
    stats.push_back(cLat);


    // max dist from starting location
    stats.push_back(
            (*std::max_element(distFromStart.begin(),distFromStart.end()))
            * fD.eDistSd);


    return stats;
}


template <class M>
void runBayesP(const std::string & file, const int obs,
        const int numSamples,const int numBurn,
        const int numStats,
        const bool save,
        const bool edgeToEdge){

    std::string edgeExt;
    if(edgeToEdge) {
        edgeExt = "edge";
    } else {
        edgeExt = "spatial";
    }

    njm::resetSeed(0);

    typedef System<M,M> S;

    S sObs("obsData.txt");
    sObs.setEdgeToEdge(edgeToEdge);

    std::vector<std::vector<int> > h;
    h = sObs.sD.history;
    h.push_back(sObs.sD.status);

    std::vector<std::string> names = {"n_inf","n_inf_2007","n_inf_2008",
                                      "n_inf_2009","n_inf_2010","n_inf_2011",
                                      "n_inf_2012","n_inf_2013","mean_year",
                                      "mean_long","mean_lat",
                                      "mean_dist_from_start",
                                      "min_long","min_lat",
                                      "max_long","max_lat",
                                      "max_dist_from_start"};

    if(obs){
        njm::toFile(names,njm::sett.datExt("obsStats_"+edgeExt+"_",".txt"),
                std::ios_base::out);
        njm::toFile(getStats(h,sObs.sD,sObs.fD),
                njm::sett.datExt("obsStats_"+edgeExt+"_",".txt"));
    }

    sObs.modelGen_r.mcmc.load(sObs.sD.history,sObs.sD.status,sObs.fD);
    sObs.modelGen_r.mcmc.sample(numSamples,numBurn,true);

    { // mean
        sObs.modelGen_r.mcmc.samples.setMean();
        std::vector<double> par = sObs.modelGen_r.mcmc.samples.getPar();
        sObs.modelGen_r.putPar(par.begin());

        if(save)
            sObs.modelGen_r.save();

        std::vector< std::vector<double> > stats;

        S s;
        s.setEdgeToEdge(edgeToEdge);
        Starts starts("startingLocations.txt");
        s.modelGen_r = sObs.modelGen_r;
        s.modelEst_r = s.modelGen_r;
        int r,t,R,T;
        R = numStats;
        T = sObs.sD.time;
        std::vector< std::vector<double> > avgInf;
        for(r = 0; r < R; ++r){
            s.modelGen_r.mcmc.samples.setRand();
            par = s.modelGen_r.mcmc.samples.getPar();
            s.modelGen_r.putPar(par.begin());
            s.modelEst_r.putPar(par.begin());

            s.reset(starts[r]);

            for(t = 0; t < T; ++t)
                s.nextPoint();
            h = s.sD.history;
            h.push_back(s.sD.status);

            // store average infections
            if (avgInf.size() == 0) {
                // init
                avgInf.resize(h.size());
                for (int i = 0; i < h.size(); ++i) {
                    for (int j = 0; j < s.fD.numNodes; ++j) {
                        if (h.at(i).at(j) >= 2) {
                            avgInf.at(i).push_back(1./R);
                        } else {
                            avgInf.at(i).push_back(0.);
                        }
                    }
                }
            } else {
                for (int i = 0; i < h.size(); ++i) {
                    for (int j = 0; j < s.fD.numNodes; ++j) {
                        if (h.at(i).at(j) >= 2) {
                            avgInf.at(i).at(j) += 1./R;
                        }
                    }
                }
            }

            stats.push_back(getStats(h,s.sD,s.fD));
        }

        njm::toFile(names,njm::sett.datExt("sampStats_mean_"+file+"_"
                        +edgeExt+"_",".txt"),std::ios_base::out);
        njm::toFile(njm::toString(stats,"\n",""),
                njm::sett.datExt("sampStats_mean_"+file+"_"+edgeExt+"_",
                        ".txt"));
        njm::toFile(njm::toString(avgInf,"\n",""),
                njm::sett.datExt("sampStats_mean_avgInf_"+file+"_"+edgeExt+"_",
                        ".txt"));
    }


    { // mle
        sObs.modelGen_r.setType(MLE);
        sObs.modelGen_r.fit(sObs.sD,sObs.tD,sObs.fD,sObs.dD,0);

        njm::toFile(njm::toString(sObs.modelGen_r.getPar()," ","\n"),
                njm::sett.datExt("sampStats_"+file+"_MLE_"+edgeExt+"_",".txt"));

        std::vector<double> par;

        std::vector< std::vector<double> > stats;

        S s;
        s.setEdgeToEdge(edgeToEdge);
        Starts starts("startingLocations.txt");
        s.modelGen_r = sObs.modelGen_r;
        s.modelGen_r.setFisher(sObs.sD,sObs.tD,sObs.fD,sObs.dD);
        s.modelEst_r = s.modelGen_r;
        int r,t,R,T;
        R = numStats;
        T = sObs.sD.time;
        std::vector< std::vector<double> > avgInf;
        for(r = 0; r < R; ++r){
            s.modelGen_r.sample(true);
            std::vector<double> par = s.modelGen_r.getPar();
            s.modelEst_r.putPar(par.begin());

            s.reset(starts[r]);

            for(t = 0; t < T; ++t)
                s.nextPoint();
            h = s.sD.history;
            h.push_back(s.sD.status);


            // store average infections
            if (avgInf.size() == 0) {
                // init
                avgInf.resize(h.size());
                for (int i = 0; i < h.size(); ++i) {
                    for (int j = 0; j < s.fD.numNodes; ++j) {
                        if (h.at(i).at(j) >= 2) {
                            avgInf.at(i).push_back(1./R);
                        } else {
                            avgInf.at(i).push_back(0.);
                        }
                    }
                }
            } else {
                for (int i = 0; i < h.size(); ++i) {
                    for (int j = 0; j < s.fD.numNodes; ++j) {
                        if (h.at(i).at(j) >= 2) {
                            avgInf.at(i).at(j) += 1./R;
                        }
                    }
                }
            }


            stats.push_back(getStats(h,s.sD,s.fD));
        }

        njm::toFile(names,njm::sett.datExt("sampStats_mle_"+file+"_"
                        +edgeExt+"_",".txt"),std::ios_base::out);
        njm::toFile(njm::toString(stats,"\n",""),
                njm::sett.datExt("sampStats_mle_"+file+"_"+edgeExt+"_",".txt"));
        njm::toFile(njm::toString(avgInf,"\n",""),
                njm::sett.datExt("sampStats_mle_avgInf_"+file+"_"+edgeExt+"_",
                        ".txt"));

    }



    // param estimates
    std::vector<std::vector<double> > parSamp;
    int i;
    for(i = 0; i < (numSamples-numBurn); ++i){
        sObs.modelGen_r.mcmc.samples.setPar(i);
        njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.getPar(),
                        " ","\n"),
                njm::sett.datExt("sampStats_"+file+"_param_"+edgeExt+"_",
                        ".txt"),
                std::ios_base::app);
    }

    for(i = 0; i < numBurn; ++i){
        sObs.modelGen_r.mcmc.samples.setPar(i,true);
        njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.getPar(),
                        " ","\n"),
                njm::sett.datExt("sampStats_"+file+"_paramBurn_"+edgeExt+"_",
                        ".txt"),
                std::ios_base::app);
    }

    // posterior mode
    sObs.modelGen_r.mcmc.samples.setMode();
    njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.getPar()," ","\n"),
            njm::sett.datExt("sampStats_"+file+"_paramMode_"+edgeExt+"_",
                    ".txt"));

    // posterior mean
    sObs.modelGen_r.mcmc.samples.setMean();
    njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.getPar()," ","\n"),
            njm::sett.datExt("sampStats_"+file+"_paramMean_"+edgeExt+"_",
                    ".txt"));

    // likelihood
    njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.ll,"\n",""),
            njm::sett.datExt("sampStats_"+file+"_ll_"+edgeExt+"_",".txt"));

    // likelihood at mean
    njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.llPt,"\n"),
            njm::sett.datExt("sampStats_"+file+"_llPt_"+edgeExt+"_",".txt"));

    // pD
    njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.pD,"\n"),
            njm::sett.datExt("sampStats_"+file+"_pD_"+edgeExt+"_",".txt"));

    // Dbar
    njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.Dbar,"\n"),
            njm::sett.datExt("sampStats_"+file+"_Dbar_"+edgeExt+"_",".txt"));

    // DIC
    njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.DIC,"\n"),
            njm::sett.datExt("sampStats_"+file+"_DIC_"+edgeExt+"_",".txt"));

}



template <class M>
void runBayesPOos(const std::string & file, const int obs,
        const int numSamples,const int numBurn,
        const int numStats,
        const bool edgeToEdge){
    // bayes P out of sample

    std::string edgeExt;
    if(edgeToEdge) {
        edgeExt = "edge";
    } else {
        edgeExt = "spatial";
    }

    njm::resetSeed(0);

    typedef System<M,M> S;

    S sObs("obsData_0-5.txt");
    sObs.setEdgeToEdge(edgeToEdge);

    std::vector<std::string> names = {"n_inf","n_inf_2012",
                                      "n_inf_2013","mean_year",
                                      "mean_long","mean_lat",
                                      "average_dist_from_start",
                                      "min_long","min_lat",
                                      "max_long","max_lat",
                                      "max_dist_from_start"};

    std::vector<std::vector<int> > baseHObs;
    {
        std::vector<int> baseHObs_raw;
        std::vector<int> add;
        njm::fromFile(baseHObs_raw,njm::sett.srcExt("obsData_0-5.txt"));
        std::vector<int>::const_iterator it,end;
        it = baseHObs_raw.begin();
        end = baseHObs_raw.end();

        CHECK_EQ(baseHObs_raw.size(),6*sObs.fD.numNodes);

        int count = 0;
        while(it != end) {
            add.push_back(*it);
            ++it;
            ++count;
            if(count == sObs.fD.numNodes) {
                baseHObs.push_back(add);
                add.clear();
                count = 0;
            }
        }
        CHECK_EQ(count,0);
        CHECK_EQ(baseHObs.size(),6);
    }

    std::vector< std::vector<int> > newHObs;
    {
        std::vector<int> newHObs_raw;
        std::vector<int> add;
        njm::fromFile(newHObs_raw,njm::sett.srcExt("obsData_6-7.txt"));
        std::vector<int>::const_iterator it,end;
        it = newHObs_raw.begin();
        end = newHObs_raw.end();

        int count = 0;
        while(it != end) {
            add.push_back(*it);
            ++it;
            ++count;
            if(count == sObs.fD.numNodes) {
                newHObs.push_back(add);
                add.clear();
                count = 0;
            }
        }
        CHECK_EQ(count,0);
        CHECK_EQ(newHObs.size(),2);
    }


    if(obs){
        njm::toFile(names,njm::sett.datExt("obsStats_Oos_"+edgeExt+"_",".txt"),
                std::ios_base::out);
        njm::toFile(getStatsOos(baseHObs,newHObs,sObs.fD),
                njm::sett.datExt("obsStats_Oos_"+edgeExt+"_",".txt"));
    }

    sObs.modelGen_r.mcmc.load(sObs.sD.history,sObs.sD.status,sObs.fD);
    sObs.modelGen_r.mcmc.sample(numSamples,numBurn,true);

    { // mean
        sObs.modelGen_r.mcmc.samples.setMean();
        std::vector<double> par = sObs.modelGen_r.mcmc.samples.getPar();
        sObs.modelGen_r.putPar(par.begin());

        std::vector< std::vector<double> > stats;

        S s;
        s.setEdgeToEdge(edgeToEdge);
        std::vector<int> startVec;
        for (int i = 0; i < sObs.fD.numNodes; ++i) {
            if(baseHObs.at(baseHObs.size() - 1).at(i) >= 2)
                startVec.push_back(i);
        }
        Starts starts(startVec);
        s.modelGen_r = sObs.modelGen_r;
        s.modelEst_r = s.modelGen_r;
        int r,t,R,T;
        R = numStats;
        T = newHObs.size();
        std::vector< std::vector<double> > avgInf;
        for(r = 0; r < R; ++r){
            s.modelGen_r.mcmc.samples.setRand();
            par = s.modelGen_r.mcmc.samples.getPar();
            s.modelGen_r.putPar(par.begin());
            s.modelEst_r.putPar(par.begin());

            s.reset(starts[r]);

            // DON'T include current status in newHSim.  This should
            // only include new observations.  The current status is
            // included in baseHObs.
            std::vector<std::vector<int> > newHSim;

            for(t = 0; t < T; ++t){
                s.nextPoint();
                newHSim.push_back(s.sD.status);
            }

            // store average infections
            if (avgInf.size() == 0) {
                // init
                avgInf.resize(newHSim.size());
                for (int i = 0; i < newHSim.size(); ++i) {
                    for (int j = 0; j < s.fD.numNodes; ++j) {
                        if (newHSim.at(i).at(j) >= 2) {
                            avgInf.at(i).push_back(1./R);
                        } else {
                            avgInf.at(i).push_back(0.);
                        }
                    }
                }
            } else {
                for (int i = 0; i < newHSim.size(); ++i) {
                    for (int j = 0; j < s.fD.numNodes; ++j) {
                        if (newHSim.at(i).at(j) >= 2) {
                            avgInf.at(i).at(j) += 1./R;
                        }
                    }
                }
            }


            CHECK_EQ(newHSim.size(),2);
            stats.push_back(getStatsOos(baseHObs,newHSim,s.fD));
        }

        njm::toFile(names,njm::sett.datExt("sampStats_mean_Oos_"+file+
                        "_"+edgeExt+"_",".txt"),std::ios_base::out);
        njm::toFile(njm::toString(stats,"\n",""),
                njm::sett.datExt("sampStats_mean_Oos_"+file+"_"+edgeExt+"_",
                        ".txt"));

        njm::toFile(njm::toString(avgInf,"\n",""),
                njm::sett.datExt("sampStats_mean_Oos_avgInf_"+file+"_"+
                        edgeExt+"_",".txt"));

    }


    { // mle
        sObs.modelGen_r.setType(MLE);
        sObs.modelGen_r.fit(sObs.sD,sObs.tD,sObs.fD,sObs.dD,0);

        njm::toFile(njm::toString(sObs.modelGen_r.getPar()," ","\n"),
                njm::sett.datExt("sampStats_Oos_"+file+"_MLE_"+edgeExt+"_",
                        ".txt"));

        std::vector<double> par;

        std::vector< std::vector<double> > stats;

        S s;
        s.setEdgeToEdge(edgeToEdge);
        std::vector<int> startVec;
        for (int i = 0; i < sObs.fD.numNodes; ++i) {
            if(baseHObs.at(baseHObs.size() - 1).at(i) >= 2)
                startVec.push_back(i);
        }
        Starts starts(startVec);
        s.modelGen_r = sObs.modelGen_r;
        s.modelGen_r.setFisher(sObs.sD,sObs.tD,sObs.fD,sObs.dD);
        s.modelEst_r = s.modelGen_r;
        int r,t,R,T;
        R = numStats;
        T = newHObs.size();
        std::vector< std::vector<double> > avgInf;
        for(r = 0; r < R; ++r){
            s.modelGen_r.sample(true);
            std::vector<double> par = s.modelGen_r.getPar();
            s.modelEst_r.putPar(par.begin());

            s.reset(starts[r]);

            // DON'T include current status in newHSim.  This should
            // only include new observations.  The current status is
            // included in baseHObs.
            std::vector<std::vector<int> > newHSim;

            for(t = 0; t < T; ++t){
                s.nextPoint();
                newHSim.push_back(s.sD.status);
            }

            // store average infections
            if (avgInf.size() == 0) {
                // init
                avgInf.resize(newHSim.size());
                for (int i = 0; i < newHSim.size(); ++i) {
                    for (int j = 0; j < s.fD.numNodes; ++j) {
                        if (newHSim.at(i).at(j) >= 2) {
                            avgInf.at(i).push_back(1./R);
                        } else {
                            avgInf.at(i).push_back(0.);
                        }
                    }
                }
            } else {
                for (int i = 0; i < newHSim.size(); ++i) {
                    for (int j = 0; j < s.fD.numNodes; ++j) {
                        if (newHSim.at(i).at(j) >= 2) {
                            avgInf.at(i).at(j) += 1./R;
                        }
                    }
                }
            }


            CHECK_EQ(newHSim.size(),2);
            stats.push_back(getStatsOos(baseHObs,newHSim,s.fD));
        }

        njm::toFile(names,njm::sett.datExt("sampStats_mle_Oos_"+file+"_"
                        +edgeExt+"_",".txt"),std::ios_base::out);
        njm::toFile(njm::toString(stats,"\n",""),
                njm::sett.datExt("sampStats_mle_Oos_"+file+"_"+edgeExt+"_",
                        ".txt"));

        njm::toFile(njm::toString(avgInf,"\n",""),
                njm::sett.datExt("sampStats_mle_Oos_avgInf_"+file+"_"
                        +edgeExt+"_",".txt"));

    }



    // param estimates
    std::vector<std::vector<double> > parSamp;
    int i;
    for(i = 0; i < (numSamples-numBurn); ++i){
        sObs.modelGen_r.mcmc.samples.setPar(i);
        njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.getPar(),
                        " ","\n"),
                njm::sett.datExt("sampStats_Oos_"+file+"_param_"+edgeExt+"_",
                        ".txt"),
                std::ios_base::app);
    }

    for(i = 0; i < numBurn; ++i){
        sObs.modelGen_r.mcmc.samples.setPar(i,true);
        njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.getPar(),
                        " ","\n"),
                njm::sett.datExt("sampStats_Oos_"+file+"_paramBurn_"+edgeExt+"_"
                        ,".txt"),
                std::ios_base::app);
    }

    // posterior mode
    sObs.modelGen_r.mcmc.samples.setMode();
    njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.getPar()," ","\n"),
            njm::sett.datExt("sampStats_Oos_"+file+"_paramMode_"+edgeExt+"_",
                    ".txt"));

    // posterior mean
    sObs.modelGen_r.mcmc.samples.setMean();
    njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.getPar()," ","\n"),
            njm::sett.datExt("sampStats_Oos_"+file+"_paramMean_"+edgeExt+"_",
                    ".txt"));

    // likelihood
    njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.ll,"\n",""),
            njm::sett.datExt("sampStats_Oos_"+file+"_ll_"+edgeExt+"_",
                    ".txt"));

    // likelihood at mean
    njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.llPt,"\n"),
            njm::sett.datExt("sampStats_Oos_"+file+"_llPt_"+edgeExt+"_",
                    ".txt"));

    // pD
    njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.pD,"\n"),
            njm::sett.datExt("sampStats_Oos_"+file+"_pD_"+edgeExt+"_",".txt"));

    // Dbar
    njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.Dbar,"\n"),
            njm::sett.datExt("sampStats_Oos_"+file+"_Dbar_"+edgeExt+"_",
                    ".txt"));

    // DIC
    njm::toFile(njm::toString(sObs.modelGen_r.mcmc.samples.DIC,"\n"),
            njm::sett.datExt("sampStats_Oos_"+file+"_DIC_"+edgeExt+"_",".txt"));

}



int main(int argc, char ** argv){
    InitGoogleLogging(argv[0]);
    ParseCommandLineFlags(&argc,&argv,true);
    if(!FLAGS_dryRun) {
        njm::sett.setup(std::string(argv[0]),FLAGS_srcDir);


        int numSamples = 100000, numBurn = 20000, numStats = 50000;
        // int numSamples = 50000, numBurn = 25000, numStats = 10000;
        // int numSamples = 20000, numBurn = 10000, numStats = 10000;
        // int numSamples = 100, numBurn = 50, numStats = 50;
        // int numSamples = 10, numBurn = 5, numStats = 5;

        if(FLAGS_save){
            std::cout << "Setup to SAVE parameters." << std::endl;
        }
        else{
            std::cout << "Setup to NOT SAVE parameters." << std::endl;
        }

#pragma omp parallel sections                   \
    shared(numSamples,numBurn,numStats)
        {
// #pragma omp section
//       {
//         runBayesP<ModelGravityGDist
//                   >("gravity",1,
//                     numSamples,numBurn,numStats,save,false);
//       }

#pragma omp section
            {
                runBayesP<Model2GravityEDist
                          >("gravity2",1,
                                  numSamples,numBurn,numStats,FLAGS_save,false);
            }

#pragma omp section
            {
                runBayesP<Model2EdgeToEdge
                          >("edgeToEdge2",0,
                                  numSamples,numBurn,numStats,FLAGS_save,true);
            }

#pragma omp section
            {
                runBayesPOos<Model2GravityEDist
                             >("gravity2",1,
                                     numSamples,numBurn,numStats,false);
            }

#pragma omp section
            {
                runBayesPOos<Model2EdgeToEdge
                             >("edgeToEdge2",0,
                                     numSamples,numBurn,numStats,true);
            }

// #pragma omp section
//     {
//       runBayesP<ModelGravityGDistTrend
// 		>("gravityTrend",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelGravityGDistTrendPow
// 		>("gravityTrendPow",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelGravityGDistTrendPowCon
// 		>("gravityTrendPowCon",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeGDist
// 		>("timeInf",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeGDistTrend
// 		>("timeInfTrend",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeGDistTrendPow
// 		>("timeInfTrendPow",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeGDistTrendPowCon
// 		>("timeInfTrendPowCon",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeExpGDist
// 		>("timeInfExp",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeExpGDistTrend
// 		>("timeInfExpTrend",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeExpGDistTrendPow
// 		>("timeInfExpTrendPow",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeExpGDistTrendPowCon
// 		>("timeInfExpTrendPowCon",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeExpCavesGDist
// 		>("timeInfExpCaves",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeExpCavesGDistTrend
// 		>("timeInfExpCavesTrend",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeExpCavesGDistTrendPow
// 		>("timeInfExpCavesTrendPow",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelTimeExpCavesGDistTrendPowCon
// 		>("timeInfExpCavesTrendPowCon",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelRad
// 		>("rad",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelGDist
// 		>("gDist",0,
// 		  numSamples,numBurn,numStats);
//     }

// #pragma omp section
//     {
//       runBayesP<ModelGDistPow
// 		>("gDistPow",0,
// 		  numSamples,numBurn,numStats);
//     }

        }

        // njm::sett.clean();

    }
    return 0;
}
