#include "test.hpp"
#include <omp.h>


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  {
    typedef ModelTimeExpCavesGDistTrendPowCon MG;

    typedef MG ME;

    typedef System<MG,ME> S;

    typedef NoTrt<ME> NT;
    typedef ProximalGDistAgent<ME> PA;
    typedef MyopicAgent<ME> MA;

    typedef WnsFeatures3<ME> F;
    typedef RankAgent<F,ME> RA;

    typedef M1SpOptim<S,RA,ME> SPO;

    typedef VanillaRunner<S,NT> R_NT;
    typedef VanillaRunner<S,PA> R_PA;
    typedef FitOnlyRunner<S,MA> R_MA;
    typedef OptimRunner<S,RA,SPO> R_RA;


    // S s;
    S s("obsData.txt");
    s.modelGen_r.setType(MLES);
    s.modelEst_r.setType(MLES);

    int numReps = 96;
    Starts starts("startingLocations.txt");

    NT nt;
    PA pa;
    MA ma;
    RA ra;

    ra.tp.jitterScale = -1;

    SPO spo;
    spo.tp.fixSample = 1;

    R_NT r_nt;
    R_PA r_pa;
    R_MA r_ma;
    R_RA r_ra;


    RunStats rs;

    njm::message(std::string("Gravity Model with time since infected ") +
		 std::string("and time dependent intercept"));

    // 15
    njm::message("15 years");

    s.fD.forecastFlat = false;
    rs = r_nt.run(s,nt,numReps,s.fD.trtStart+15,starts);
    njm::message("  No Trt: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

    rs = r_pa.run(s,pa,numReps,s.fD.trtStart+15,starts);
    njm::message("Proximal: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

    rs = r_ma.run(s,ma,numReps,s.fD.trtStart+15,starts);
    njm::message("  Myopic: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

    // 15 flat
    njm::message("15 years flat");

    s.fD.forecastFlat = true;
    rs = r_nt.run(s,nt,numReps,s.fD.trtStart+15,starts);
    njm::message("  No Trt: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

    rs = r_pa.run(s,pa,numReps,s.fD.trtStart+15,starts);
    njm::message("Proximal: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

    rs = r_ma.run(s,ma,numReps,s.fD.trtStart+15,starts);
    njm::message("  Myopic: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

    // 25
    njm::message("25 years");

    s.fD.forecastFlat = false;
    rs = r_nt.run(s,nt,numReps,s.fD.trtStart+25,starts);
    njm::message("  No Trt: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

    rs = r_pa.run(s,pa,numReps,s.fD.trtStart+25,starts);
    njm::message("Proximal: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

    rs = r_ma.run(s,ma,numReps,s.fD.trtStart+25,starts);
    njm::message("  Myopic: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

    // 25 flat
    njm::message("25 years flat");

    s.fD.forecastFlat = true;
    rs = r_nt.run(s,nt,numReps,s.fD.trtStart+25,starts);
    njm::message("  No Trt: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

    rs = r_pa.run(s,pa,numReps,s.fD.trtStart+25,starts);
    njm::message("Proximal: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

    rs = r_ma.run(s,ma,numReps,s.fD.trtStart+25,starts);
    njm::message("  Myopic: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

  }


  {
    typedef ModelTimeExpCavesGDist MG;

    typedef MG ME;

    typedef System<MG,ME> S;

    typedef NoTrt<ME> NT;
    typedef ProximalGDistAgent<ME> PA;
    typedef MyopicAgent<ME> MA;

    typedef WnsFeatures3<ME> F;
    typedef RankAgent<F,ME> RA;

    typedef M1SpOptim<S,RA,ME> SPO;

    typedef VanillaRunner<S,NT> R_NT;
    typedef VanillaRunner<S,PA> R_PA;
    typedef FitOnlyRunner<S,MA> R_MA;
    typedef OptimRunner<S,RA,SPO> R_RA;


    // S s;
    S s("obsData.txt");
    s.modelGen_r.setType(MLES);
    s.modelEst_r.setType(MLES);

    int numReps = 96;
    Starts starts("startingLocations.txt");

    NT nt;
    PA pa;
    MA ma;
    RA ra;

    ra.tp.jitterScale = -1;

    SPO spo;
    spo.tp.fixSample = 1;

    R_NT r_nt;
    R_PA r_pa;
    R_MA r_ma;
    R_RA r_ra;


    RunStats rs;

    njm::message("Gravity Model with time since infected");

    // 15
    njm::message("15 years");

    s.fD.forecastFlat = false;
    rs = r_nt.run(s,nt,numReps,s.fD.trtStart+15,starts);
    njm::message("  No Trt: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

    rs = r_pa.run(s,pa,numReps,s.fD.trtStart+15,starts);
    njm::message("Proximal: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

    rs = r_ma.run(s,ma,numReps,s.fD.trtStart+15,starts);
    njm::message("  Myopic: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

    // 25
    njm::message("25 years");

    s.fD.forecastFlat = false;
    rs = r_nt.run(s,nt,numReps,s.fD.trtStart+25,starts);
    njm::message("  No Trt: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

    rs = r_pa.run(s,pa,numReps,s.fD.trtStart+25,starts);
    njm::message("Proximal: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

    rs = r_ma.run(s,ma,numReps,s.fD.trtStart+25,starts);
    njm::message("  Myopic: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

  }


  {
    typedef ModelGravityGDist MG;

    typedef MG ME;

    typedef System<MG,ME> S;

    typedef NoTrt<ME> NT;
    typedef ProximalGDistAgent<ME> PA;
    typedef MyopicAgent<ME> MA;

    typedef WnsFeatures3<ME> F;
    typedef RankAgent<F,ME> RA;

    typedef M1SpOptim<S,RA,ME> SPO;

    typedef VanillaRunner<S,NT> R_NT;
    typedef VanillaRunner<S,PA> R_PA;
    typedef FitOnlyRunner<S,MA> R_MA;
    typedef OptimRunner<S,RA,SPO> R_RA;


    // S s;
    S s("obsData.txt");
    s.modelGen_r.setType(MLES);
    s.modelEst_r.setType(MLES);

    int numReps = 96;
    Starts starts("startingLocations.txt");

    NT nt;
    PA pa;
    MA ma;
    RA ra;

    ra.tp.jitterScale = -1;

    SPO spo;
    spo.tp.fixSample = 1;

    R_NT r_nt;
    R_PA r_pa;
    R_MA r_ma;
    R_RA r_ra;


    RunStats rs;


    njm::message("Gravity Model");

    // 15
    njm::message("15 years");

    s.fD.forecastFlat = false;
    rs = r_nt.run(s,nt,numReps,s.fD.trtStart+15,starts);
    njm::message("  No Trt: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

    rs = r_pa.run(s,pa,numReps,s.fD.trtStart+15,starts);
    njm::message("Proximal: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

    rs = r_ma.run(s,ma,numReps,s.fD.trtStart+15,starts);
    njm::message("  Myopic: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

    // 25
    njm::message("25 years");

    s.fD.forecastFlat = false;
    rs = r_nt.run(s,nt,numReps,s.fD.trtStart+25,starts);
    njm::message("  No Trt: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

    rs = r_pa.run(s,pa,numReps,s.fD.trtStart+25,starts);
    njm::message("Proximal: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

    rs = r_ma.run(s,ma,numReps,s.fD.trtStart+25,starts);
    njm::message("  Myopic: "
		 + njm::toString(rs.smean(),"")
		 + "  (" + njm::toString(rs.seMean(),"") + ")");

  }

  njm::sett.clean();
  return 0;
}
