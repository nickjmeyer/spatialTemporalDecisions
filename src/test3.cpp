#include "test3.hpp"

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);

  typedef GravityModel GM;
  typedef GravityParam GP;
  
  typedef System<GM,GP,GM,GP> S;
  
  typedef NoTrt<GM,GP> NA;
  typedef ProximalAgent<GM,GP> PA;

  typedef PlainRunner<S,NA> R_NA;
  typedef PlainRunner<S,PA> R_PA;

  S s;
  NA na;
  PA pa;

  R_NA r_na;
  R_PA r_pa;

  njm::message(r_na.run(s,na,1000,15));
  njm::message(r_pa.run(s,pa,1000,15));

  njm::sett.clean();
  return 0;
}
