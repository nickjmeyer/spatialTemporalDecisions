#include "sgdTune.hpp"


Design::Design(){
  numVar=0;
}


void Design::add(const std::string var, const double val, const double inc,
		 const double lower = std::numeric_limits<double>::min(),
		 const double upper = std::numeric_limits<double>::max()){
  numVar++;
  vars.push_back(var);
  levs[var] = std::vector<double>();
  int i,line=1;
  double cand;
  for(i = -line; i <= line; i++)
    if(i != 0){
      cand = val + ((double)i)*inc;
      if(cand >= upper)
	levs[var].push_back(upper - .001);
      else if(cand <= lower)
	levs[var].push_back(lower + .001);
      else
	levs[var].push_back(cand);
    }
  levs[var].erase(std::unique(levs[var].begin(),levs[var].end()),
		  levs[var].end());
}


int Design::next(){
  std::string var;
  int i,ind;
  for(i=0; i<numVar; i++){
    var = vars.at(i);
    ind = ++inds[var];
    if(ind == (int)levs[var].size())
      ind = 0;
    inds[var] = ind;
    vals[var] = levs[var].at(ind);
    if(ind > 0)
      return 1;
  }
  return 0;
}


void Design::reset(){
  int i;
  for(i=0; i<numVar; i++){
    inds[vars.at(i)] = 0;
    vals[vars.at(i)] = levs[vars.at(i)].at(0);
  }
}


void Design::clear(){
  numVar = 0;
  vars.clear();
  inds.clear();
  vals.clear();
  levs.clear();
}


double Design::operator[](const std::string var) const{
  return vals.at(var);
}


int main(int argc, char ** argv){
  njm::sett.set(argc,argv);
  
#pragma omp parallel
  {
    System<GravityModel,GravityParam> s;
    s.estParam_r = s.genParam_r;
    s.reset();

    RankToyAgent<GravityModel,GravityParam> rA;
  
    PlainRunner<System,RankToyAgent,GravityModel,GravityParam> pR;
    M1SgdOptim<System,RankToyAgent,GravityModel,GravityParam> sO;


    Design d;
    double jitter,rateDecay,rateStart,cur,can;
    int numYears = 11;

    jitter = .2;
    rateDecay = .975;
    rateStart = 10;
    
    sO.tp.numYears = numYears;

    rA.tp.weights.ones();
    sO.optim(s,rA);
    cur=pR.run(s,rA,100,numYears);

    double scale=1.0;
    int r,i,R=250,change;
    arma::mat res(R,4 + rA.numFeatures);
    std::vector<double> weights = rA.tp.getPar();

    res(0,0) = cur;
    res(0,1) = jitter;
    res(0,2) = rateDecay;
    res(0,3) = rateStart;
    for(i=0; i<rA.numFeatures; i++)
      res(0,4+i) = weights.at(i);    
    
    for(r=1; r<R; r++){
      d.clear();
      d.add("jitter", jitter, 0.02*scale, 0.0, 2.0);
      d.add("rateDecay", rateDecay, 0.005*scale, 0.5, 1.0);
      d.add("rateStart", rateStart, 2.0*scale, 0.0, 100.0);
      d.reset();

      scale*=.99;

      change=0;
      do{
	sO.tp.jitter = d["jitter"];
	sO.tp.rateDecay = d["rateDecay"];
	sO.tp.rate = d["rateStart"];
	
	rA.tp.weights.ones();
	sO.optim(s,rA);	
	can=pR.run(s,rA,100,numYears);

	if(can < cur){
	  change=1;
	  jitter = d["jitter"];
	  rateDecay = d["rateDecay"];
	  rateStart = d["rateStart"];
	  weights = rA.tp.getPar();
	  cur = can;
	}
      } while(d.next());

      if(!change){
	sO.tp.jitter = jitter;
	sO.tp.rateDecay = rateDecay;
	sO.tp.rate = rateStart;
	
	rA.tp.weights.ones();
	sO.optim(s,rA);		
	cur=pR.run(s,rA,100,numYears);
      }
      
      res(r,0) = cur;
      res(r,1) = jitter;
      res(r,2) = rateDecay;
      res(r,3) = rateStart;
      for(i=0; i<rA.numFeatures; i++)
	res(r,4+i) = weights.at(i);

#pragma omp critical
      {      
	res.save(njm::sett.datExt("res"+njm::toString(omp_get_thread_num(),
						      "",0,0),".txt"),
		 arma::raw_ascii);
      }
      
      if(omp_get_thread_num()==0)
	std::cout << "Thread 0 completed " << njm::toString(r,"",4,0)
		  << " out of " << R << "\r" << std::flush;
    }
    if(omp_get_thread_num()==0)
      std::cout << std::endl;
  }


  njm::sett.clean();
  return 0;
}
