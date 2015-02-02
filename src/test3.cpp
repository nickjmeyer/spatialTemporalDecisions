#include "test3.hpp"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);
  
  System<RangeModel,RangeParam,
	 RangeModel,RangeParam> s;

  NoTrt<RangeModel,RangeParam> nT;

  PlainRunner<System<RangeModel,RangeParam,
		     RangeModel,RangeParam>,
	      NoTrt<RangeModel,RangeParam> > pR;

  njm::message(pR.run(s,nT,1000,15));
	      

  
  njm::sett.clean();
  return 0;
}
