#include "mcmcRadius.hpp"


enum parInd{INTCP_=0,RADIUS_=1,TRTP_=2,TRTA_=3};

void RadiusSamples::setMean(){
  intcpSet = radiusSet = trtPreSet = trtActSet = 0.0;

  intcpSet = std::accumulate(intcp.begin(),intcp.end(),0.0);
  intcpSet /= double(numSamples);
  radiusSet = std::accumulate(radius.begin(),radius.end(),0.0);
  radiusSet /= double(numSamples);
  trtPreSet = std::accumulate(trtPre.begin(),trtPre.end(),0.0);
  trtPreSet /= double(numSamples);
  trtActSet = std::accumulate(trtAct.begin(),trtAct.end(),0.0);
  trtActSet /= double(numSamples);
}


void RadiusSamples::setRand(){
  intcpSet = radiusSet = trtPreSet = trtActSet = 0.0;

  int i = njm::runifInterv(0,numSamples);

  intcpSet = intcp.at(i);
  radiusSet = radius.at(i);
  trtPreSet = trtPre.at(i);
  trtActSet = trtAct.at(i);

}


std::vector<double> RadiusSamples::getPar() const {
  std::vector<double> par;
  par.push_back(intcpSet);
  par.push_back(radiusSet);
  par.push_back(trtActSet);
  par.push_back(trtPreSet);

  return par;
}



void RadiusMcmc::load(const std::vector<std::vector<int> > & history,
		      const std::vector<int> & status,
		      const FixedData & fD){
  std::vector<std::vector<int> > all;
  all = history;
  all.push_back(status);
  load(all,fD);
}



void RadiusMcmc::load(const std::vector<std::vector<int> > & history,
		      const FixedData & fD){
  numNodes=fD.numNodes;
  T=(int)history.size();
  numCovar=fD.numCovar;
  samples.numCovar = numCovar;

  priorTrtMean = fD.priorTrtMean;

  infHist.resize(numNodes*T);
  trtPreHist.resize(numNodes*T);
  trtActHist.resize(numNodes*T);
  d = fD.dist;
  cc.resize(numNodes*numNodes);
  covar = fD.covar;
  timeInf.resize(numNodes*T);
  int i,j;
  for(i = 0; i < numNodes; ++i){
    for(j = 0; j < T; ++j){// get the histories of infection and treatments
      infHist.at(i*T + j)=(history.at(j).at(i) < 2 ? 0 : 1);
      trtPreHist.at(i*T + j)=(history.at(j).at(i) == 1 ? 1 : 0);
      trtActHist.at(i*T + j)=(history.at(j).at(i) == 3 ? 1 : 0);
    }

    // while you're looping, get d and cc
    for(j = 0; j < numNodes; ++j){
      cc.at(i*numNodes + j)=fD.caves.at(i)*fD.caves.at(j);
    }
  }


  int val;
  for(i = 0; i < numNodes; ++i){
    val = 0;
    for(j = 0; j < T; ++j){
      if(infHist.at(i*T + j) == 1)
	++val;
      timeInf.at(i*T + j) = val;
    }
  }

}

void RadiusMcmc::sample(int const numSamples, int const numBurn,
			const bool saveBurn){
  std::vector<double> par = {-3.0, // intcp
			     100, // radius
			     0.0, // trtAct
			     0.0}; // trtPre
  sample(numSamples,numBurn,par,saveBurn);
}


void RadiusMcmc::sample(int const numSamples, int const numBurn,
			const std::vector<double> & par,
			const bool saveBurn){
  samples.numSamples = numSamples - numBurn;
  samples.numBurn = numBurn;

  // priors
  int thin=1;
  double intcp_mean=0,intcp_var=100,
    radius_mean=100,radius_var=100,
    trtPre_mean=priorTrtMean,trtPre_var=1,
    trtAct_mean=priorTrtMean,trtAct_var=1;


  int i,j;
  // set containers for current and candidate samples
  std::vector<double>::const_iterator it = par.begin();
  intcp_cur=intcp_can= *it++;
  radius_cur=radius_can= *it++;
  trtPre_cur=trtPre_can= *it++;
  trtAct_cur=trtAct_can= *it++;

  // set containers for storing all non-burned samples
  samples.intcp.clear();
  samples.intcpHist.clear();
  samples.intcp.reserve(numSamples-numBurn);
  samples.intcpHist.reserve(numBurn);

  samples.radius.clear();
  samples.radiusHist.clear();
  samples.radius.reserve(numSamples-numBurn);
  samples.radiusHist.reserve(numBurn);

  samples.trtPre.clear();
  samples.trtPreHist.clear();
  samples.trtPre.reserve(numSamples-numBurn);
  samples.trtPreHist.reserve(numBurn);

  samples.trtAct.clear();
  samples.trtActHist.clear();
  samples.trtAct.reserve(numSamples-numBurn);
  samples.trtActHist.reserve(numBurn);


  samples.ll.clear();
  samples.llHist.clear();
  samples.ll.reserve(numSamples-numBurn);
  samples.llHist.reserve(numBurn);


  // get the likelihood with the current parameters
  ll_cur=ll_can=ll();

  // set the MH tuning parameters
  acc=att= std::vector<int>(5,0);
  mh=std::vector<double>(5,0.5);
  // tau=std::vector<double>(numCovar+2,0.0);

  // mu=std::vector<double>(numCovar+2,0.0);
  // mu.at(numCovar+INTCP_) = -3;

  double upd;
  double R;

  int displayOn=1;
  int display=0;

  // do a bunch of nonsense...
  for(i=0; i<numSamples; ++i){
    if(display && i%displayOn==0){
      printf("SLM...%6s: %6d\r","iter",i);
      fflush(stdout);
    }


    // sample intcp
    ++att.at(INTCP_);
    upd=intcp_cur+mh.at(INTCP_)*njm::rnorm01();
    intcp_can=upd;

    // get new likelihood
    ll_can=ll();


    R=ll_can + (-.5/intcp_var)*std::pow(intcp_can - intcp_mean,2.0)
      - ll_cur - (-.5/intcp_var)*std::pow(intcp_cur - intcp_mean,2.0);


    // accept?
    if(std::log(njm::runif01()) < R){
      ++acc.at(INTCP_);
      intcp_cur=intcp_can;
      ll_cur=ll_can;
    }
    else{
      intcp_can=intcp_cur;
      ll_can=ll_cur;
    }


    // sample trtPre
    ++att.at(TRTP_);
    upd=trtPre_cur+mh.at(TRTP_)*njm::rnorm01();
    trtPre_can=upd;


    // get new likelihood
    ll_can=ll();

    R=ll_can + (-.5/trtPre_var)*std::pow(trtPre_can - trtPre_mean,2.0)
      - ll_cur - (-.5/trtPre_var)*std::pow(trtPre_cur - trtPre_mean,2.0);

    // accept?
    if(std::log(njm::runif01()) < R){
      ++acc.at(TRTP_);
      trtPre_cur=trtPre_can;
      ll_cur=ll_can;
    }
    else{
      trtPre_can=trtPre_cur;
      ll_can=ll_cur;
    }



    // sample trtAct
    ++att.at(TRTA_);
    upd=trtAct_cur+mh.at(TRTA_)*njm::rnorm01();
    trtAct_can=upd;

    // get new likelihood
    ll_can=ll();

    R=ll_can + (-.5/trtAct_var)*std::pow(trtAct_can - trtAct_mean,2.0)
      - ll_cur - (-.5/trtAct_var)*std::pow(trtAct_cur - trtAct_mean,2.0);


    // accept?
    if(std::log(njm::runif01()) < R){
      ++acc.at(TRTA_);
      trtAct_cur=trtAct_can;
      ll_cur=ll_can;
    }
    else{
      trtAct_can=trtAct_cur;
      ll_can=ll_cur;
    }




    // sample radius
    ++att.at(RADIUS_);
    upd=std::exp(std::log(radius_cur)+mh.at(RADIUS_)*njm::rnorm01());
    radius_can=upd;


    // get new likelihood
    ll_can=ll();


    R=ll_can + (-.5/radius_var)*std::pow(std::log(radius_can)
					 - radius_mean,2.0)
      - ll_cur - (-.5/radius_var)*std::pow(std::log(radius_cur)
					   - radius_mean,2.0);


    // accept?
    if(std::log(njm::runif01()) < R){
      ++acc.at(RADIUS_);
      radius_cur=radius_can;
      ll_cur=ll_can;
    }
    else{
      radius_can=radius_cur;
      ll_can=ll_cur;
    }





    if(i<numBurn){
      // time for tuning!
      int len=int(mh.size());
      double accRatio;
      for(j = 0; j < len; ++j){
	if(att.at(j) > 50){
	  accRatio=((double)acc.at(j))/((double)att.at(j));
	  if(accRatio < .3)
	    mh.at(j)*=.8;
	  else if(accRatio > .6)
	    mh.at(j)*=1.2;

	  acc.at(j)=0;
	  att.at(j)=0;
	}
      }
      if(saveBurn){
	samples.intcpHist.push_back(intcp_cur);
	samples.radiusHist.push_back(radius_cur);
	samples.trtPreHist.push_back(trtPre_cur);
	samples.trtActHist.push_back(trtAct_cur);
      }
    }
    else if(i%thin==0){
      // save the samples
      samples.intcp.push_back(intcp_cur);
      samples.radius.push_back(radius_cur);
      samples.trtPre.push_back(trtPre_cur);
      samples.trtAct.push_back(trtAct_cur);

      samples.ll.push_back(ll_cur);
    }
  }

  // get likelihood evaluated at posterior mean
  samples.setMean();
  intcp_can = samples.intcpSet;
  radius_can = samples.radiusSet;
  trtPre_can = samples.trtPreSet;
  trtAct_can = samples.trtActSet;

  samples.llPt = ll();

  samples.Dbar=-2.0*std::accumulate(samples.ll.begin(),
				    samples.ll.end(),
				    0.0)/double(numSamples);
  samples.pD=samples.Dbar - -2.0*samples.llPt;
  samples.DIC=samples.pD + samples.Dbar;


  if(display)
    printf("\33[2K\r");
}



double RadiusMcmc::ll(){
  int i,j,k;
  double llVal,wontProb,prob,expProb,baseProb,baseProbInit;

  llVal = 0.0;
  for(i=1; i<T; i++){// loop over time interval that has changed
    for(j=0; j<numNodes; j++){
      if(infHist.at(j*T + i-1)==0){// if county is susceptible get infProb
	wontProb=1.0;
	// set a base number to decrease floating point operations
	if(trtPreHist.at(j*T + i-1)==0)
	  baseProbInit=intcp_can;
	else
	  baseProbInit=intcp_can - trtPre_can;

	for(k=0; k<numNodes; k++){
	  // if county is infected it affects the infProb
	  if(infHist.at(k*T + i-1)==1){
	    // calculate infProb
	    baseProb=baseProbInit;
	    if(d.at(j*numNodes + k) > radius_can)
	      baseProb -= 500.0;

	    if(trtActHist.at(k*T + i-1)==1)
	      baseProb -= trtAct_can;

	    expProb=std::exp(baseProb);

	    wontProb*=1.0/(1.0+expProb);
	  }
	}

	prob=1.0-wontProb;

	if(!(prob > 0.0))
	  prob=std::exp(-30.0);
	else if(!(prob < 1.0))
	  prob=1.0 - std::exp(-30.0);

	if(infHist.at(j*T + i)==0)
	  llVal+=std::log(1-prob);
	else
	  llVal+=std::log(prob);
      }
    }
  }

  return llVal;
}
