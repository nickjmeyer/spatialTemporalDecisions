#include "test.hpp"


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef Model2GPowGDist M;
  typedef System<M,M> S;


  S s;
  s.modelGen_r.setType(MLES);
  s.modelEst_r.setType(MLES);

  std::vector<std::string> names = {"power"};
  double power = s.modelGen_r.getPar(names)[0];
  s.modelGen_r.setPar("power",std::log(power));

  std::vector<double> par = s.modelGen_r.getPar();
  s.modelEst_r.putPar(par.begin());
  s.revert();

  Starts starts(15,s.fD.numNodes);
  s.reset(starts[3]);


  int i;
  for(i = 0; i < s.fD.trtStart; ++i){
    s.nextPoint();
  }

  std::cout << "numInfected: " << s.sD.numInfected << std::endl;

  s.modelEst.fit(s.sD,s.tD,s.fD,s.dD,false);

  std::ofstream ofs;
  ofs.open("hessTestData.txt",std::ofstream::out);
  ofs << "par0, par1, val0, val1, llVal, llGrad, llHess\n";

  std::vector<double> changeVals;
  int numChange = 20;
  int pi,pj;
  std::vector<double> origPar = s.modelGen.getPar();
  std::vector<double> newPar;
  for(pi = 0; pi < int(s.modelGen.numPars); ++pi){
    for(pj = 0; pj < int(s.modelGen.numPars); ++pj){
      std::cout << std::setfill(' ') << std::setw(2) << pi
		<< ", " << std::setfill(' ') << std::setw(2) << pj
		<< "\r" << std::flush;

      for(i = 0; i < (numChange+1); ++i){
	double change = 0.01 * 2.0 *(double(i)/double(numChange) - 0.5);
	newPar = origPar;
	newPar.at(pj) += change;
	s.modelGen.putPar(newPar.begin());

	double llVal = s.modelGen.logll(s.sD,s.tD,s.fD,s.dD);
	std::vector<double> llGrad = s.modelGen.logllGrad(s.sD,s.tD,s.fD,s.dD);
	std::vector<double> llHess = s.modelGen.logllHess(s.sD,s.tD,s.fD,s.dD);

	ofs << pi << ", "
	    << pj << ", "
	    << newPar.at(pi) << ", "
	    << newPar.at(pj) << ", "
	    << llVal << ", "
	    << llGrad[pi] << ", "
	    << llHess[pi*s.modelGen.numPars + pj] << "\n";
      }

    }
  }
  std::cout << std::endl;

  ofs.close();

  // std::cout << "gen: " << njm::toString(s.modelGen.getPar()) << std::endl;
  // std::cout << "est: " << njm::toString(s.modelEst.getPar()) << std::endl;

  // for(i = 0; i < 10; i++){
  //   s.modelEst.sample();
  //   std::cout << i << ": " << njm::toString(s.modelEst.getPar()) << std::endl;
  //   std::cout << "==========================================" << std::endl;
  // }

  njm::sett.clean();
  return 0;
}
