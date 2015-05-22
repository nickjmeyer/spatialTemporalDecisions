#include "runM1Mles.hpp"


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef ModelTimeExpCavesGDistTrendPowCon MG;
  
  typedef MG ME;

  typedef System<MG,ME> S;

  typedef NoTrt<ME> NT;
  typedef ProximalGDistAgent<ME> PA;
  typedef MyopicAgent<ME> MA;
  
  typedef ToyFeatures5<ME> F;
  typedef RankAgent<F,ME> RA;
  typedef OsspAgent<ME> OA;

  typedef M1SpOptim<S,RA,ME> SPO;
  typedef M1OsspOptim<S,OA,F,ME> OSSPO;

  typedef VanillaRunner<S,NT> R_NT;
  typedef VanillaRunner<S,PA> R_PA;
  typedef FitOnlyRunner<S,MA> R_MA;
  typedef OptimRunner<S,RA,SPO> R_RA;
  typedef OptimRunner<S,OA,OSSPO> R_OA;


  S s;
  s.modelGen_r.setType(MLE);
  s.modelEst_r.setType(MLE);

  int numReps = 96;
  Starts starts(numReps,s.fD.numNodes);

  NT nt;
  PA pa;
  MA ma;
  RA ra;
  OA oa;

  SPO spo;
  OSSPO osspo;

  R_NT r_nt;
  R_PA r_pa;
  R_MA r_ma;
  R_RA r_ra;
  R_OA r_oa;


  RunStats rs;

  // njm::message("par: " + njm::toString(s.modelGen_r.getPar()," ",""));


  // std::vector<double> d = s.fD.gDist;

  // std::vector<double> gamma = {1,2,4};
  // std::vector<double> K = {3,5,10};

  // std::vector<double>::iterator gIt,gEnd;
  // std::vector<double>::iterator kIt,kEnd;

  // for(gIt = gamma.begin(); gIt != gEnd; ++gIt){
  //   for(kIt = K.begin(); kIt != kEnd; ++kIt){

  // s.fD.gDist = d;
  // std::for_each(s.fD.gDist.begin(),s.fD.gDist.end(),
  // 		[&gIt,&kIt](double & x){x = std::pow(x,*gIt)/ *kIt;});

  // double gamma,k;
  // njm::fromFile(gamma,njm::sett.srcExt("gamma.txt"));
  // njm::fromFile(k,njm::sett.srcExt("k.txt"));

  // s.fD.gDist = d;
  // std::for_each(s.fD.gDist.begin(),s.fD.gDist.end(),
  // 		[&gamma,&k](double & x){
  // 		  x = std::pow(x/k,gamma);
  // 		});

  // s.preCompData();

  // s.modelGen_r = MG(s.fD);
  // s.modelGen_r.read();
  // s.modelEst_r = ME(s.fD);
  // s.modelGen_r.setType(MLES);
  // s.modelEst_r.setType(MLES);

  // std::cout << std::endl << std::endl
  // 	    << "k: " << *kIt << std::endl
  // 	    << "gamma: " << *gIt << std::endl;
  
  // std::cout << std::endl << std::endl
  // 	    << "k: " << k << std::endl
  // 	    << "gamma: " << gamma << std::endl;
  

  rs = r_nt.run(s,nt,numReps,s.fD.finalT,starts);
  njm::message("   No treatment: "
	       + njm::toString(rs.smean(),"")
	       + "  (" + njm::toString(rs.seMean(),"") + ")");
  
  rs = r_pa.run(s,pa,numReps,s.fD.finalT,starts);
  njm::message("       Proximal: "
  	       + njm::toString(rs.smean(),"")
  	       + "  (" + njm::toString(rs.seMean(),"") + ")");
  
  rs = r_ma.run(s,ma,numReps,s.fD.finalT,starts);
  njm::message("         Myopic: "
  	       + njm::toString(rs.smean(),"")
  	       + "  (" + njm::toString(rs.seMean(),"") + ")");
  
  rs = r_ra.run(s,ra,spo,numReps,s.fD.finalT,starts);
  njm::message("  Policy Search: "
  	       + njm::toString(rs.smean(),"")
  	       + "  (" + njm::toString(rs.seMean(),"") + ")");

  // osspo.tp.N = 100;
  // osspo.tp.jitterScale = 4.0;
  // osspo.name = "M1Ossp_" + njm::toString(osspo.tp.N,"",0,0)
  // 	+ "_" + njm::toString(osspo.tp.jitterScale,"",0,0);
  rs = r_oa.run(s,oa,osspo,numReps,s.fD.finalT,starts);
  njm::message("One Step Polish: "
  		   + njm::toString(rs.smean(),"")
  		   + "  (" + njm::toString(rs.seMean(),"") + ")");
  
  // osspo.tp.N = 100;
  // osspo.tp.jitterScale = 1.0;
  // osspo.name = "M1Ossp_" + njm::toString(osspo.tp.N,"",0,0)
  // 	+ "_" + njm::toString(osspo.tp.jitterScale,"",0,0);
  // rs = r_oa.run(s,oa,osspo,numReps,s.fD.finalT,starts);
  // njm::message("One Step Polish: "
  // 		   + njm::toString(rs.smean(),"")
  // 		   + "  (" + njm::toString(rs.seMean(),"") + ")");
  
  // osspo.tp.N = 1000;
  // osspo.tp.jitterScale = 4.0;
  // osspo.name = "M1Ossp_" + njm::toString(osspo.tp.N,"",0,0)
  // 	+ "_" + njm::toString(osspo.tp.jitterScale,"",0,0);
  // rs = r_oa.run(s,oa,osspo,numReps,s.fD.finalT,starts);
  // njm::message("One Step Polish: "
  // 		   + njm::toString(rs.smean(),"")
  // 		   + "  (" + njm::toString(rs.seMean(),"") + ")");
  
  // osspo.tp.N = 1000;
  // osspo.tp.jitterScale = 1.0;
  // osspo.name = "M1Ossp_" + njm::toString(osspo.tp.N,"",0,0)
  // 	+ "_" + njm::toString(osspo.tp.jitterScale,"",0,0);
  // rs = r_oa.run(s,oa,osspo,numReps,s.fD.finalT,starts);
  // njm::message("One Step Polish: "
  // 		   + njm::toString(rs.smean(),"")
  // 		   + "  (" + njm::toString(rs.seMean(),"") + ")");
  
  //   }
  // }
  
  return 0;
}

