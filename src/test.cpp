#include "test.hpp"
#include "omp.h"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  {
    typedef GravityModel GM;
    typedef GravityParam GP;
    typedef GravityModel EM;
    typedef GravityParam EP;

    typedef System<GM,GP,EM,EP> S;
    typedef NoTrt<EM,EP> NT;

    typedef PlainRunner<S,NT> PR;
  
    S s;

    NT nt;
    PR pr;

    std::vector<double> scales{
      0.001,0.005,0.01,0.05,0.1,0.5,
	1.0,1.5,2.0,2.5,3.0,3.5,
	4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5};


    std::vector<double> parOrig = s.paramGen_r.getPar();
    
    int i;
    for(i = 0; i < (int)scales.size(); i++){
      std::vector<double> par = parOrig;
      double scale = scales.at(i);
      std::for_each(par.begin(),par.end(),
      		    [&scale](double & x){x*=scale;});
      s.paramGen_r.putPar(par);
      
      s.paramEst_r = s.paramGen_r;
      s.reset();

      printf("[% 8.6f] --> (% 6.4f)\n",scale,pr.run(s,nt,300,15));
    }
  }


  njm::sett.clean();
  return 0;
}
