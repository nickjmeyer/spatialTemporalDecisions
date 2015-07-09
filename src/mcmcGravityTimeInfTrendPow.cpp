#include "mcmcGravityTimeInfTrendPow.hpp"


enum parInd{INTCP_=0,ALPHA_=1,POWER_=2,TREND_=3,TRENDPOW_=4,XI_=5,
	    TRTP_=6,TRTA_=7};

void GravityTimeInfTrendPowSamples::setMean(){
  intcpSet = alphaSet = powerSet = trendSet = trendPowSet = xiSet =
    trtPreSet = trtActSet = 0.0;
  betaSet.resize(numCovar);
  std::fill(betaSet.begin(),betaSet.end(),0.0);

  intcpSet = std::accumulate(intcp.begin(),intcp.end(),0.0);
  intcpSet /= double(numSamples);
  alphaSet = std::accumulate(alpha.begin(),alpha.end(),0.0);
  alphaSet /= double(numSamples);
  powerSet = std::accumulate(power.begin(),power.end(),0.0);
  powerSet /= double(numSamples);
  trendSet = std::accumulate(trend.begin(),trend.end(),0.0);
  trendSet /= double(numSamples);
  trendPowSet = std::accumulate(trendPow.begin(),trendPow.end(),0.0);
  trendPowSet /= double(numSamples);
  xiSet = std::accumulate(xi.begin(),xi.end(),0.0);
  xiSet /= double(numSamples);
  trtPreSet = std::accumulate(trtPre.begin(),trtPre.end(),0.0);
  trtPreSet /= double(numSamples);
  trtActSet = std::accumulate(trtAct.begin(),trtAct.end(),0.0);
  trtActSet /= double(numSamples);

  int j = 0;
  std::for_each(beta.begin(),beta.end(),
		[this,&j](const double & x){
		  betaSet.at(j++ % numCovar) += x;
		});
  std::for_each(betaSet.begin(),betaSet.end(),
		[this](double & x){x /= double(numSamples);});
}


void GravityTimeInfTrendPowSamples::setRand(){
  intcpSet = alphaSet = powerSet = trendSet = trendPowSet = xiSet =
    trtPreSet = trtActSet = 0.0;
  betaSet.resize(numCovar);
  std::fill(betaSet.begin(),betaSet.end(),0.0);

  int i = njm::runifInterv(0,numSamples);

  intcpSet = intcp.at(i);
  alphaSet = alpha.at(i);
  powerSet = power.at(i);
  trendSet = trend.at(i);
  trendPowSet = trendPow.at(i);
  xiSet = xi.at(i);
  trtPreSet = trtPre.at(i);
  trtActSet = trtAct.at(i);

  int j = 0;
  std::for_each(betaSet.begin(),betaSet.end(),
		[this,&i,&j](double & x){
		  x = beta.at(i*numCovar + j++);});
}


void GravityTimeInfTrendPowSamples::setPar(const int i,const bool fromBurn){
  intcpSet = alphaSet = powerSet = trendSet = trendPowSet = xiSet =
    trtPreSet = trtActSet = 0.0;
  betaSet.resize(numCovar);
  std::fill(betaSet.begin(),betaSet.end(),0.0);

  if(fromBurn){
    intcpSet = intcpHist.at(i);
    alphaSet = alphaHist.at(i);
    powerSet = powerHist.at(i);
    trendSet = trendHist.at(i);
    trendPowSet = trendPowHist.at(i);
    xiSet = xiHist.at(i);
    trtPreSet = trtPreHist.at(i);
    trtActSet = trtActHist.at(i);

    int j = 0;
    std::for_each(betaSet.begin(),betaSet.end(),
		  [this,&i,&j](double & x){
		    x = betaHist.at(i*numCovar + j++);});
  }
  else{
    intcpSet = intcp.at(i);
    alphaSet = alpha.at(i);
    powerSet = power.at(i);
    trendSet = trend.at(i);
    trendPowSet = trendPow.at(i);
    xiSet = xi.at(i);
    trtPreSet = trtPre.at(i);
    trtActSet = trtAct.at(i);

    int j = 0;
    std::for_each(betaSet.begin(),betaSet.end(),
		  [this,&i,&j](double & x){
		    x = beta.at(i*numCovar + j++);});
  }
}


std::vector<double> GravityTimeInfTrendPowSamples::getPar() const {
  std::vector<double> par = {intcpSet};
  par.insert(par.end(),betaSet.begin(),betaSet.end());
  par.push_back(alphaSet);
  par.push_back(powerSet);
  par.push_back(trendSet);
  par.push_back(trendPowSet);
  par.push_back(xiSet);
  par.push_back(trtActSet);
  par.push_back(trtPreSet);

  return par;
}



void GravityTimeInfTrendPowMcmc
::load(const std::vector<std::vector<int> > & history,
       const std::vector<int> & status,
       const FixedData & fD){
  std::vector<std::vector<int> > all;
  all = history;
  all.push_back(status);
  load(all,fD);
}



void GravityTimeInfTrendPowMcmc
::load(const std::vector<std::vector<int> > & history,
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
  cc.resize(numNodes*numNodes);
  covar = fD.covar;
  timeInfMinOne.resize(numNodes*T);
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
      timeInfMinOne.at(i*T + j) = val - 1.0;
    }
  }

}


void GravityTimeInfTrendPowMcmc::sample(int const numSamples,
					int const numBurn,
					const bool saveBurn){
  std::vector<double> beta (numCovar,0.0);
  std::vector<double> par = {-3.0, // intcp
			     0.1, // alpha
			     0.1, // power
			     0.1, // trend
			     0.1, // trendPow
			     0.1, // xi
			     0.0, // trtAct
			     0.0}; // trtPre
  par.insert(par.begin()+1,beta.begin(),beta.end());
  sample(numSamples,numBurn,par,saveBurn);
}



void GravityTimeInfTrendPowMcmc::sample(int const numSamples, int const numBurn,
					const std::vector<double> & par,
					const bool saveBurn){
  samples.numSamples = numSamples - numBurn;
  samples.numBurn = numBurn;

  // priors
  int thin=1;
  double intcp_mean=0,intcp_var=100,beta_mean=0,beta_var=10,alpha_mean=0,
    alpha_var=1,power_mean=0,power_var=1,xi_mean=0,xi_var=1.0,
    trend_mean=0,trend_var=1,
    trendPow_mean=0,trendPow_var=1,
    trtPre_mean=priorTrtMean,trtPre_var=1,
    trtAct_mean=priorTrtMean,trtAct_var=1;


  int i,j;
  // set containers for current and candidate samples
  std::vector<double>::const_iterator it = par.begin();
  intcp_cur=intcp_can= *it++;

  beta_cur.clear();
  for(i = 0; i < numCovar; ++i)
    beta_cur.push_back(*it++);
  beta_can = beta_cur;

  alpha_cur=alpha_can= (*it < 0.00001 ? 0.01 : *it);
  ++it;
  power_cur=power_can= (*it < 0.00001 ? 0.01 : *it);
  ++it;
  trend_cur=trend_can=*it++;
  trendPow_cur=trendPow_can=*it++;
  xi_cur=xi_can= (*it < 0.00001 ? 0.01 : *it);
  ++it;
  trtPre_cur=trtPre_can= *it++;
  trtAct_cur=trtAct_can= *it++;

  // set containers for storing all non-burned samples
  samples.intcp.clear();
  samples.intcpHist.clear();
  samples.intcp.reserve(numSamples-numBurn);
  samples.intcpHist.reserve(numBurn);

  samples.beta.clear();
  samples.betaHist.clear();
  samples.beta.reserve((numSamples-numBurn)*numCovar);
  samples.betaHist.reserve((numBurn)*numCovar);

  samples.alpha.clear();
  samples.alphaHist.clear();
  samples.alpha.reserve(numSamples-numBurn);
  samples.alphaHist.reserve(numBurn);

  samples.power.clear();
  samples.powerHist.clear();
  samples.power.reserve(numSamples-numBurn);
  samples.powerHist.reserve(numBurn);

  samples.trend.clear();
  samples.trendHist.clear();
  samples.trend.reserve(numSamples-numBurn);
  samples.trendHist.reserve(numBurn);

  samples.trendPow.clear();
  samples.trendPowHist.clear();
  samples.trendPow.reserve(numSamples-numBurn);
  samples.trendPowHist.reserve(numBurn);

  samples.xi.clear();
  samples.xiHist.clear();
  samples.xi.reserve(numSamples-numBurn);
  samples.xiHist.reserve(numBurn);

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


  covarBeta_cur.resize(numNodes);
  updateCovarBeta(covarBeta_cur,covar,beta_cur,numNodes,numCovar);
  covarBeta_can = covarBeta_cur;

  alphaW_cur.resize(numNodes*numNodes);
  updateAlphaW(alphaW_cur,d,cc,alpha_cur,power_cur,numNodes);
  alphaW_can = alphaW_cur;

  xiTimeInfMinOne_cur = timeInfMinOne;
  updateXiTimeInf(xiTimeInfMinOne_cur,xi_cur);
  xiTimeInfMinOne_can = xiTimeInfMinOne_cur;

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
      printf("McmcGravityTimeInfTrendPow...%6s: %6d\r","iter",i);
      fflush(stdout);
    }


    // sample intcp
    ++att.at(numCovar+INTCP_);
    upd=intcp_cur+mh.at(numCovar+INTCP_)*njm::rnorm01();
    intcp_can=upd;

    // get new likelihood
    ll_can=ll();


    R=ll_can + (-.5/intcp_var)*std::pow(intcp_can - intcp_mean,2.0)
      - ll_cur - (-.5/intcp_var)*std::pow(intcp_cur - intcp_mean,2.0);


    // accept?
    if(std::log(njm::runif01()) < R){
      ++acc.at(numCovar+INTCP_);
      intcp_cur=intcp_can;
      ll_cur=ll_can;
    }
    else{
      intcp_can=intcp_cur;
      ll_can=ll_cur;
    }



    // sample beta
    for(j = 0; j < numCovar; ++j){
      ++att.at(j);

      upd=beta_cur.at(j)+mh.at(j)*njm::rnorm01();
      beta_can.at(j)=upd;

      updateCovarBeta(covarBeta_can,covar,
		      beta_cur.at(j),beta_can.at(j),
		      j,numCovar);

      // get new likelihood
      ll_can=ll();

      R=ll_can + (-.5/beta_var)*std::pow(beta_can.at(j) - beta_mean,2.0)
	- ll_cur - (-.5/beta_var)*std::pow(beta_cur.at(j) - beta_mean,2.0);

      // accept?
      if(std::log(njm::runif01()) < R){
	++acc.at(j);
	beta_cur.at(j)=beta_can.at(j);
	ll_cur=ll_can;
	covarBeta_cur=covarBeta_can;
      }
      else{
	beta_can.at(j)=beta_cur.at(j);
	covarBeta_can=covarBeta_cur;
	ll_can=ll_cur;
      }
    }


    // sample trtPre
    ++att.at(numCovar+TRTP_);
    upd=trtPre_cur+mh.at(numCovar+TRTP_)*njm::rnorm01();
    trtPre_can=upd;


    // get new likelihood
    ll_can=ll();

    R=ll_can + (-.5/trtPre_var)*std::pow(trtPre_can - trtPre_mean,2.0)
      - ll_cur - (-.5/trtPre_var)*std::pow(trtPre_cur - trtPre_mean,2.0);

    // accept?
    if(std::log(njm::runif01()) < R){
      ++acc.at(numCovar+TRTP_);
      trtPre_cur=trtPre_can;
      ll_cur=ll_can;
    }
    else{
      trtPre_can=trtPre_cur;
      ll_can=ll_cur;
    }



    // sample trtAct
    ++att.at(numCovar+TRTA_);
    upd=trtAct_cur+mh.at(numCovar+TRTA_)*njm::rnorm01();
    trtAct_can=upd;

    // get new likelihood
    ll_can=ll();

    R=ll_can + (-.5/trtAct_var)*std::pow(trtAct_can - trtAct_mean,2.0)
      - ll_cur - (-.5/trtAct_var)*std::pow(trtAct_cur - trtAct_mean,2.0);


    // accept?
    if(std::log(njm::runif01()) < R){
      ++acc.at(numCovar+TRTA_);
      trtAct_cur=trtAct_can;
      ll_cur=ll_can;
    }
    else{
      trtAct_can=trtAct_cur;
      ll_can=ll_cur;
    }




    // sample alpha
    ++att.at(numCovar+ALPHA_);

    logAlpha_cur=std::log(alpha_cur);

    upd=std::exp(logAlpha_cur + mh.at(numCovar+ALPHA_)*njm::rnorm01());
    alpha_can=upd;
    logAlpha_can=std::log(alpha_can);

    // update alphaW
    updateAlphaW(alphaW_can,alpha_cur,alpha_can,numNodes);

    // get new likelihood
    ll_can=ll();


    R=ll_can + (-.5/alpha_var)*std::pow(logAlpha_can - alpha_mean,2.0)
      - ll_cur - (-.5/alpha_var)*std::pow(logAlpha_cur - alpha_mean,2.0);

    // accept?
    if(std::log(njm::runif01()) < R){
      ++acc.at(numCovar+ALPHA_);
      alpha_cur=alpha_can;
      alphaW_cur=alphaW_can;
      ll_cur=ll_can;
    }
    else{
      alpha_can=alpha_cur;
      alphaW_can=alphaW_cur;
      ll_can=ll_cur;
    }




    // sample power
    ++att.at(numCovar+POWER_);
    upd=std::exp(std::log(power_cur)+mh.at(numCovar+POWER_)*njm::rnorm01());
    power_can=upd;

    // update alphaW
    updateAlphaW(alphaW_can,d,cc,alpha_cur,power_can,numNodes);

    // get new likelihood
    ll_can=ll();


    R=ll_can + (-.5/power_var)*std::pow(std::log(power_can) - power_mean,2.0)
      - ll_cur - (-.5/power_var)*std::pow(std::log(power_cur) - power_mean,2.0);


    // accept?
    if(std::log(njm::runif01()) < R){
      ++acc.at(numCovar+POWER_);
      power_cur=power_can;
      alphaW_cur=alphaW_can;
      ll_cur=ll_can;
    }
    else{
      power_can=power_cur;
      alphaW_can=alphaW_cur;
      ll_can=ll_cur;
    }



    // sample trend
    ++att.at(numCovar+TREND_);
    upd=trend_cur+mh.at(numCovar+TREND_)*njm::rnorm01();
    trend_can=upd;

    // get new likelihood
    ll_can=ll();


    R=ll_can + (-.5/trend_var)*std::pow(trend_can - trend_mean,2.0)
      - ll_cur - (-.5/trend_var)*std::pow(trend_cur - trend_mean,2.0);


    // accept?
    if(std::log(njm::runif01()) < R){
      ++acc.at(numCovar+TREND_);
      trend_cur=trend_can;
      ll_cur=ll_can;
    }
    else{
      trend_can=trend_cur;
      ll_can=ll_cur;
    }




    // sample trendPow
    ++att.at(numCovar+TRENDPOW_);
    upd=trendPow_cur+mh.at(numCovar+TRENDPOW_)*njm::rnorm01();
    trendPow_can=upd;

    // get new likelihood
    ll_can=ll();


    R=ll_can + (-.5/trendPow_var)*std::pow(trendPow_can - trendPow_mean,2.0)
      - ll_cur - (-.5/trendPow_var)*std::pow(trendPow_cur - trendPow_mean,2.0);


    // accept?
    if(std::log(njm::runif01()) < R){
      ++acc.at(numCovar+TRENDPOW_);
      trendPow_cur=trendPow_can;
      ll_cur=ll_can;
    }
    else{
      trendPow_can=trendPow_cur;
      ll_can=ll_cur;
    }




    // sample xi
    ++att.at(numCovar+XI_);
    upd=xi_cur+mh.at(numCovar+XI_)*njm::rnorm01();
    xi_can=upd;

    updateXiTimeInf(xiTimeInfMinOne_can,xi_can/xi_cur);

    // get new likelihood
    ll_can=ll();


    R=ll_can + (-.5/xi_var)*std::pow(xi_can - xi_mean,2.0)
      - ll_cur - (-.5/xi_var)*std::pow(xi_cur - xi_mean,2.0);



    // accept?
    if(std::log(njm::runif01()) < R){
      ++acc.at(numCovar+XI_);
      xi_cur=xi_can;
      xiTimeInfMinOne_cur = xiTimeInfMinOne_can;
      ll_cur=ll_can;
    }
    else{
      xi_can=xi_cur;
      xiTimeInfMinOne_can=xiTimeInfMinOne_cur;
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
	samples.betaHist.insert(samples.beta.end(),
				beta_cur.begin(),
				beta_cur.end());
	samples.alphaHist.push_back(alpha_cur);
	samples.powerHist.push_back(power_cur);
	samples.trendHist.push_back(trend_cur);
	samples.trendPowHist.push_back(trendPow_cur);
	samples.xiHist.push_back(xi_cur);
	samples.trtPreHist.push_back(trtPre_cur);
	samples.trtActHist.push_back(trtAct_cur);
      }
    }
    else if(i%thin==0){
      // save the samples
      samples.intcp.push_back(intcp_cur);
      samples.beta.insert(samples.beta.end(),beta_cur.begin(),beta_cur.end());
      samples.alpha.push_back(alpha_cur);
      samples.power.push_back(power_cur);
      samples.trend.push_back(trend_cur);
      samples.trendPow.push_back(trendPow_cur);
      samples.xi.push_back(xi_cur);
      samples.trtPre.push_back(trtPre_cur);
      samples.trtAct.push_back(trtAct_cur);

      samples.ll.push_back(ll_cur);
    }
  }

  // get likelihood evaluated at posterior mean
  samples.setMean();
  intcp_can = samples.intcpSet;
  beta_can = samples.betaSet;
  alpha_can = samples.alphaSet;
  power_can = samples.powerSet;
  trend_can = samples.trendSet;
  trendPow_can = samples.trendPowSet;
  xi_can = samples.xiSet;
  trtPre_can = samples.trtPreSet;
  trtAct_can = samples.trtActSet;

  updateCovarBeta(covarBeta_can,covar,beta_can,numNodes,numCovar);
  updateAlphaW(alphaW_can,d,cc,alpha_can,power_can,numNodes);

  xiTimeInfMinOne_can = timeInfMinOne;
  updateXiTimeInf(xiTimeInfMinOne_can,xi_can);

  samples.llPt = ll();

  samples.Dbar=-2.0*std::accumulate(samples.ll.begin(),
				    samples.ll.end(),
				    0.0)/double(samples.numSamples);
  samples.pD=samples.Dbar - -2.0*samples.llPt;
  samples.DIC=samples.pD + samples.Dbar;


  if(display)
    printf("\33[2K\r");
}



double GravityTimeInfTrendPowMcmc::ll(){
  int i,j,k;
  double llVal,wontProb,prob,expProb,baseProb,baseProbInit;

  llVal = 0.0;
  for(i=1; i<T; i++){// loop over time interval that has changed
    for(j=0; j<numNodes; j++){
      if(infHist.at(j*T + i-1)==0){// if county is susceptible get infProb
	wontProb=1.0;
	// set a base number to decrease floating point operations
	if(trtPreHist.at(j*T + i-1)==0)
	  baseProbInit=intcp_can+covarBeta_can.at(j);
	else
	  baseProbInit=intcp_can+covarBeta_can.at(j) - trtPre_can;

	for(k=0; k<numNodes; k++){
	  // if county is infected it affects the infProb
	  if(infHist.at(k*T + i-1)==1){
	    // calculate infProb
	    baseProb=baseProbInit;
	    if(j < k)
	      baseProb -= alphaW_can.at(j*numNodes + k);
	    else
	      baseProb -= alphaW_can.at(k*numNodes + j);

	    // i is time + 1
	    baseProb += trend_can*std::pow(double(i),trendPow_can);

	    baseProb += xiTimeInfMinOne_can.at(k*T + i-1);

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
