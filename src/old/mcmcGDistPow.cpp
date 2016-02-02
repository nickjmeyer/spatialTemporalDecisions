#include "mcmcGDistPow.hpp"


enum parInd{INTCP_=0,ALPHA_=1,POWER_=2,TRTP_=3,TRTA_=4};



void GDistPowSamples::setMode(){
  intcpSet = DensityEst(intcp).max().second;
  alphaSet = DensityEst(alpha).max().second;
  powerSet = DensityEst(power).max().second;
  trtPreSet = DensityEst(trtPre).max().second;
  trtActSet = DensityEst(trtAct).max().second;
}



void GDistPowSamples::setMean(){
  intcpSet = alphaSet = powerSet = trtPreSet = trtActSet = 0.0;

  intcpSet = std::accumulate(intcp.begin(),intcp.end(),0.0);
  intcpSet /= double(numSamples);
  alphaSet = std::accumulate(alpha.begin(),alpha.end(),0.0);
  alphaSet /= double(numSamples);
  powerSet = std::accumulate(power.begin(),power.end(),0.0);
  powerSet /= double(numSamples);
  trtPreSet = std::accumulate(trtPre.begin(),trtPre.end(),0.0);
  trtPreSet /= double(numSamples);
  trtActSet = std::accumulate(trtAct.begin(),trtAct.end(),0.0);
  trtActSet /= double(numSamples);
}


void GDistPowSamples::setRand(){
  intcpSet = alphaSet = powerSet = trtPreSet = trtActSet = 0.0;

  int i = njm::runifInterv(0,numSamples);

  intcpSet = intcp.at(i);
  alphaSet = alpha.at(i);
  powerSet = power.at(i);
  trtPreSet = trtPre.at(i);
  trtActSet = trtAct.at(i);
}

void GDistPowSamples::setPar(const int i,const bool fromBurn){
  intcpSet = alphaSet = powerSet = trtPreSet = trtActSet = 0.0;
  
if(fromBurn){
  intcpSet = intcpBurn.at(i);
  alphaSet = alphaBurn.at(i);
  powerSet = powerBurn.at(i);
  trtPreSet = trtPreBurn.at(i);
  trtActSet = trtActBurn.at(i);
}
else{
  intcpSet = intcp.at(i);
  alphaSet = alpha.at(i);
  powerSet = power.at(i);
  trtPreSet = trtPre.at(i);
  trtActSet = trtAct.at(i);
}

}




std::vector<double> GDistPowSamples::getPar() const {
  std::vector<double> par {intcpSet};
  par.push_back(alphaSet);
  par.push_back(powerSet);
  par.push_back(trtActSet);
  par.push_back(trtPreSet);
  
  return par;
}



void GDistPowMcmc::load(const std::vector<std::vector<int> > & history,
			const std::vector<int> & status,
			const FixedData & fD){
  std::vector<std::vector<int> > all;
  all = history;
  all.push_back(status);
  load(all,fD);
}



void GDistPowMcmc::load(const std::vector<std::vector<int> > & history,
			const FixedData & fD){
  numNodes=fD.numNodes;
  T=(int)history.size();
  numCovar=fD.numCovar;
  samples.numCovar = numCovar;

  priorTrtMean = fD.priorTrtMean;
  
  infHist.resize(numNodes*T);
  trtPreHist.resize(numNodes*T);
  trtActHist.resize(numNodes*T);
  d = fD.gDist;
  timeInf.resize(numNodes*T);
  int i,j;
  for(i = 0; i < numNodes; ++i){
    for(j = 0; j < T; ++j){// get the histories of infection and treatments
      infHist.at(i*T + j)=(history.at(j).at(i) < 2 ? 0 : 1);
      trtPreHist.at(i*T + j)=(history.at(j).at(i) == 1 ? 1 : 0);
      trtActHist.at(i*T + j)=(history.at(j).at(i) == 3 ? 1 : 0);
    }

  }


  int val;
  for(i = 0; i < numNodes; ++i){
    val = 0;
    for(j = 0; j < T; ++j){
      if(infHist.at(i*T + j) == 1) // this should be 1
	++val;
      timeInf.at(i*T + j) = val;
    }
  }
  
}


void GDistPowMcmc::sample(int const numSamples, int const numBurn,
const bool saveBurn){
  std::vector<double> par = {-3.0, // intcp
			     0.1, // alpha
			     0.1, // power
			     0.0, // trtAct
			     0.0}; // trtPre
  sample(numSamples,numBurn,par,saveBurn);
}


void GDistPowMcmc::sample(int const numSamples, int const numBurn,
			  const std::vector<double> & par,
const bool saveBurn){
  samples.numSamples = numSamples - numBurn;
samples.numBurn = numBurn;
  
  // priors
  int thin=1;
  double intcp_mean=0,intcp_var=100,alpha_mean=0,
    alpha_var=1,power_mean=0,power_var=1,
    trtPre_mean=priorTrtMean,trtPre_var=1,
    trtAct_mean=priorTrtMean,trtAct_var=1;


  int i,j;
  // set containers for current and candidate samples
  std::vector<double>::const_iterator it = par.begin();
  intcp_cur=intcp_can= *it++;

  alpha_cur=alpha_can= (*it < 0.00001 ? 0.01 : *it);
  ++it;
  power_cur=power_can= 0.0;
  ++it;
  trtAct_cur=trtAct_can= *it++;
  trtPre_cur=trtPre_can= *it++;

  // set containers for storing all non-burned samples
  samples.intcp.clear();
samples.intcpBurn.clear();
  samples.intcp.reserve(numSamples-numBurn);
samples.intcpBurn.reserve(numBurn);

  samples.alpha.clear();
samples.alphaBurn.clear();
  samples.alpha.reserve(numSamples-numBurn);
samples.alphaBurn.reserve(numBurn);

  samples.power.clear();
samples.powerBurn.clear();
  samples.power.reserve(numSamples-numBurn);
samples.powerBurn.reserve(numBurn);

  samples.trtPre.clear();
samples.trtPreBurn.clear();
  samples.trtPre.reserve(numSamples-numBurn);
samples.trtPreBurn.reserve(numBurn);

  samples.trtAct.clear();
samples.trtActBurn.clear();
  samples.trtAct.reserve(numSamples-numBurn);
samples.trtActBurn.reserve(numBurn);


  samples.ll.clear();
samples.llBurn.clear();
  samples.ll.reserve(numSamples-numBurn);
samples.llBurn.reserve(numBurn);


  // get the likelihood with the current parameters
  ll_cur=ll_can=ll();

  // set the MH tuning parameters
  acc=att= std::vector<int>(par.size(),0);
  mh=std::vector<double>(par.size(),0.5);
  // tau=std::vector<double>(numCovar+2,0.0);
  
  // mu=std::vector<double>(numCovar+2,0.0);
  // mu.at(numCovar+INTCP_) = -3;
  
  double upd;
  double R;
  
  double logAlpha_cur,logAlpha_can;

  int displayOn=1;
  int display=0;

  // do a bunch of nonsense...
  for(i=0; i<numSamples; ++i){
    if(display && i%displayOn==0){
      printf("McmcGDistPow...%6s: %6d\r","iter",i);
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




    // sample alpha
    ++att.at(ALPHA_);
    
    logAlpha_cur=std::log(alpha_cur);

    upd=std::exp(logAlpha_cur + mh.at(ALPHA_)*njm::rnorm01());
    alpha_can=upd;
    logAlpha_can=std::log(alpha_can);

    // get new likelihood
    ll_can=ll();
    
    
    R=ll_can + (-.5/alpha_var)*std::pow(logAlpha_can - alpha_mean,2.0)
      - ll_cur - (-.5/alpha_var)*std::pow(logAlpha_cur - alpha_mean,2.0);

    // accept?
    if(std::log(njm::runif01()) < R){
      ++acc.at(ALPHA_);
      alpha_cur=alpha_can;
      ll_cur=ll_can;
    }
    else{
      alpha_can=alpha_cur;
      ll_can=ll_cur;
    }




    // sample power
    ++att.at(POWER_);
    upd=power_cur+mh.at(POWER_)*njm::rnorm01();
    power_can=upd;

    // get new likelihood
    ll_can=ll();


    R=ll_can + (-.5/power_var)*std::pow(power_can - power_mean,2.0)
      - ll_cur - (-.5/power_var)*std::pow(power_cur - power_mean,2.0);


    // accept?
    if(std::log(njm::runif01()) < R){
      ++acc.at(POWER_);
      power_cur=power_can;
      ll_cur=ll_can;
    }
    else{
      power_can=power_cur;
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
      samples.intcpBurn.push_back(intcp_cur);
      samples.alphaBurn.push_back(alpha_cur);
      samples.powerBurn.push_back(power_cur);
      samples.trtPreBurn.push_back(trtPre_cur);
      samples.trtActBurn.push_back(trtAct_cur);
}
    }
    else if(i%thin==0){
      // save the samples
      samples.intcp.push_back(intcp_cur);
      samples.alpha.push_back(alpha_cur);
      samples.power.push_back(power_cur);
      samples.trtPre.push_back(trtPre_cur);
      samples.trtAct.push_back(trtAct_cur);
      
      samples.ll.push_back(ll_cur);
    }
  }

  // get likelihood evaluated at posterior mean
  samples.setMean();
  intcp_can = samples.intcpSet;
  alpha_can = samples.alphaSet;
  power_can = samples.powerSet;
  trtPre_can = samples.trtPreSet;
  trtAct_can = samples.trtActSet;

  samples.llPt = ll();

  samples.Dbar=-2.0*std::accumulate(samples.ll.begin(),
				    samples.ll.end(),
				    0.0)/double(samples.numSamples);
  samples.pD=samples.Dbar - -2.0*samples.llPt;
  samples.DIC=samples.pD + samples.Dbar;


  if(display)
    printf("\33[2K\r");
}



double GDistPowMcmc::ll(){
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

	    baseProb -= alpha_can * std::pow(d[j*numNodes + k],power_can);
	    
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



