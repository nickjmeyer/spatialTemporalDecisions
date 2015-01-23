#include "mcmc.hpp"


enum parInd{INTCP_=0,ALPHA_=1,POWER_=2,TRTP_=3,TRTA_=4,XI_=5};


void GravityMcmc::load(const SimData & sD, const TrtData & tD,
		       const FixedData & fD, const DynamicData & dD){
  numNodes=fD.numNodes;
  T=(int)sD.history.size()+1;
  covar=fD.numCovar;
  
  infHist.resize(numNodes,T);
  trtPreHist.resize(numNodes,T);
  trtActHist.resize(numNodes,T);
  d.resize(numNodes,numNodes);
  cc.resize(numNodes,numNodes);
  Xcov.resize(numNodes,covar);
  int i,j;
  for(i=0; i<numNodes; i++){
    for(j=0; j<T-1; j++){// get the histories of infection and treatments
      infHist(i,j)=(sD.history.at(j).at(i) < 2 ? 0 : 1);
      trtPreHist(i,j)=(sD.history.at(j).at(i) == 1 ? 1 : 0);
      trtActHist(i,j)=(sD.history.at(j).at(i) == 3 ? 1 : 0);
    }
    // grid only has the past, now get the present
    infHist(i,j)=(sD.status.at(i) < 2 ? 0 : 1);
    trtPreHist(i,j)=(sD.status.at(i) == 1 ? 1 : 0);
    trtActHist(i,j)=(sD.status.at(i) == 3 ? 1 : 0);

    // while you're looping, get d and cc
    for(j=0; j<numNodes; j++){
      d(i,j)=fD.dist.at(i*numNodes +j);
      cc(i,j)=fD.caves.at(i)*fD.caves.at(j);
    }
    
    // might as well get Xcov too...
    for(j=0; j<covar; j++)
      Xcov(i,j)=fD.covar.at(i*covar+j);
  }

  timeInf = arma::conv_to<arma::mat>::from(arma::cumsum(infHist,1));
}



void GravityMcmc::sample(int const numSamples, int const numBurn){
  // priors
  int thin=1;
  double intcp_mean=0,intcp_var=100,beta_mean=0,beta_var=10,alpha_mean=0,
    alpha_var=1,power_mean=0,power_var=1,trtPre_mean=-4,trtPre_var=100,
    trtAct_mean=-4,trtAct_var=100;


  int i,j;
  // set containers for current and candidate samples
  intcp_cur=intcp_can=-3;
  beta_cur=beta_can=arma::zeros<arma::colvec>(covar);
  alpha_cur=alpha_can=.1;
  power_cur=power_can=.1;
  trtPre_cur=trtPre_can=0;
  trtAct_cur=trtAct_can=0;

  // set containers for storing all non-burned samples
  samples.intcp.zeros(numSamples-numBurn);
  samples.beta.zeros(numSamples-numBurn,covar);
  samples.alpha.zeros(numSamples-numBurn);
  samples.power.zeros(numSamples-numBurn);
  samples.trtPre.zeros(numSamples-numBurn);
  samples.trtAct.zeros(numSamples-numBurn);

  samples.ll.zeros(numSamples-numBurn);

  XcovBeta_cur=XcovBeta_can=arma::zeros<arma::colvec>(numNodes);
  XcovBeta_can=XcovBeta_cur=Xcov*beta_can;

  alphaW_cur.resize(numNodes,numNodes);
  alphaW_can.resize(numNodes,numNodes);
  alphaW_cur=arma::as_scalar(alpha_cur)
    *(d/arma::pow(cc,power_cur));
  alphaW_can=arma::as_scalar(alpha_can)
    *(d/arma::pow(cc,power_can));

  
  // get the likelihood with the current parameters
  llVec_cur.zeros(T-1);
  llVec_can.zeros(T-1);
  ll_cur=ll_can=ll();
  llVec_cur=llVec_can;

  // set the MH tuning parameters
  acc=att=arma::zeros<arma::colvec> (covar+5)+.5;
  mh=arma::zeros<arma::colvec> (covar+5)+.5;
  tau=arma::ones<arma::colvec> (covar+2);
  
  mu=arma::zeros<arma::colvec> (covar+2);
  mu(covar+INTCP_)=-3;
  
  double upd;
  double R;
  
  double logAlpha_cur,logAlpha_can;

  int displayOn=50;
  int display=0;

  // do a bunch of nonsense...
  for(i=0; i<numSamples; i++){
    if(display && i%displayOn==0){
      printf("SLM...%6s: %6d\r","iter",i);
      fflush(stdout);
    }


    // sample intcp
    att(covar+INTCP_)++;
    upd=intcp_cur+mh(covar+INTCP_)*njm::rnorm01();
    intcp_can=upd;
    
    // get new likelihood
    ll_can=ll();
    
    
    R=ll_can + (-.5/intcp_var)*std::pow(intcp_can - intcp_mean,2.0)
      - ll_cur - (-.5/intcp_var)*std::pow(intcp_cur - intcp_mean,2.0);
      
    
    // accept?
    if(std::log(njm::runif01()) < R){
      acc(covar+INTCP_)++;
      intcp_cur=intcp_can;
      ll_cur=ll_can;
      llVec_cur=llVec_can;
    }
    else{
      intcp_can=intcp_cur;
      ll_can=ll_cur;
      llVec_can=llVec_cur;
    }

    // sample beta
    for(j=0; j<covar; j++){
      att(j)++;
      
      upd=beta_cur(j)+mh(j)*njm::rnorm01();
      beta_can(j)=upd;
      
      XcovBeta_can=Xcov*beta_can;
      // get new likelihood
      ll_can=ll();

      R=ll_can + (-.5/beta_var)*std::pow(beta_can(j) - beta_mean,2.0)
	- ll_cur - (-.5/beta_var)*std::pow(beta_cur(j) - beta_mean,2.0);
      
      // accept?
      if(std::log(njm::runif01()) < R){
	acc(j)++;
	beta_cur(j)=beta_can(j);
	ll_cur=ll_can;
	llVec_cur=llVec_can;
	XcovBeta_cur=XcovBeta_can;
      }
      else{
	beta_can(j)=beta_cur(j);
	XcovBeta_can=XcovBeta_cur;
	ll_can=ll_cur;
	llVec_can=llVec_cur;
      }
    }


    // sample trtPre
    att(covar+TRTP_)++;
    upd=trtPre_cur+mh(covar+TRTP_)*njm::rnorm01();
    trtPre_can=upd;


    // get new likelihood
    ll_can=ll();

    R=ll_can + (-.5/trtPre_var)*std::pow(trtPre_can - trtPre_mean,2.0)
      - ll_cur - (-.5/trtPre_var)*std::pow(trtPre_cur - trtPre_mean,2.0);

    // accept?
    if(std::log(njm::runif01()) < R){
      acc(covar+TRTP_)++;
      trtPre_cur=trtPre_can;
      ll_cur=ll_can;
      llVec_cur=llVec_can;
    }
    else{
      trtPre_can=trtPre_cur;
      ll_can=ll_cur;
      llVec_can=llVec_cur;
    }



    // sample trtAct
    att(covar+TRTA_)++;
    upd=trtAct_cur+mh(covar+TRTA_)*njm::rnorm01();
    trtAct_can=upd;

    // get new likelihood
    ll_can=ll();

    R=ll_can + (-.5/trtAct_var)*std::pow(trtAct_can - trtAct_mean,2.0)
      - ll_cur - (-.5/trtAct_var)*std::pow(trtAct_cur - trtAct_mean,2.0);

    
    // accept?
    if(std::log(njm::runif01()) < R){
      acc(covar+TRTA_)++;
      trtAct_cur=trtAct_can;
      ll_cur=ll_can;
      llVec_cur=llVec_can;
    }
    else{
      trtAct_can=trtAct_cur;
      ll_can=ll_cur;
      llVec_can=llVec_cur;
    }




    // sample alpha
    att(covar+ALPHA_)++;
    
    logAlpha_cur=std::log(alpha_cur);

    upd=std::exp(logAlpha_cur + mh(covar+ALPHA_)*njm::rnorm01());
    alpha_can=upd;
    logAlpha_can=std::log(alpha_can);

    // update alphaW
    alphaW_can*=arma::as_scalar(alpha_can/alpha_cur);

    // get new likelihood
    ll_can=ll();
    
    
    R=ll_can + (-.5/alpha_var)*std::pow(logAlpha_can - alpha_mean,2.0)
      - ll_cur - (-.5/alpha_var)*std::pow(logAlpha_cur - alpha_mean,2.0);

    // accept?
    if(std::log(njm::runif01()) < R){
      acc(covar+ALPHA_)++;
      alpha_cur=alpha_can;
      alphaW_cur=alphaW_can;
      ll_cur=ll_can;
      llVec_cur=llVec_can;
    }
    else{
      alpha_can=alpha_cur;
      alphaW_can=alphaW_cur;
      ll_can=ll_cur;
      llVec_can=llVec_cur;
    }




    // sample power
    att(covar+POWER_)++;
    upd=std::exp(std::log(power_cur)+mh(covar+POWER_)*njm::rnorm01());
    power_can=upd;

    // update alphaW
    alphaW_can=alpha_can*d/arma::pow(cc,power_can);

    // get new likelihood
    ll_can=ll();


    R=ll_can + (-.5/power_var)*std::pow(std::log(power_can) - power_mean,2.0)
      - ll_cur - (-.5/power_var)*std::pow(std::log(power_cur) - power_mean,2.0);


    // accept?
    if(std::log(njm::runif01()) < R){
      acc(covar+POWER_)++;
      power_cur=power_can;
      alphaW_cur=alphaW_can;
      ll_cur=ll_can;
      llVec_cur=llVec_can;
    }
    else{
      power_can=power_cur;
      alphaW_can=alphaW_cur;
      ll_can=ll_cur;
      llVec_can=llVec_cur;
    }

    if(i<numBurn){
      // time for tuning!
      int len=mh.n_elem;
      double accRatio;
      for(j=0; j<len; j++){
	if(att(j)>50){
	  accRatio=((double)acc(j))/((double)att(j));
	  if(accRatio < .3)
	    mh(j)*=.8;
	  else if(accRatio > .6)
	    mh(j)*=1.2;
	  
	  acc(j)=0;
	  att(j)=0;
	}
      }      
    }
    else if(i%thin==0){
      // save the samples
      samples.intcp(i-numBurn)=intcp_cur;
      for(j=0; j<covar; j++)
	samples.beta(i-numBurn,j)=beta_cur(j);
      samples.alpha(i-numBurn)=alpha_cur;
      samples.power(i-numBurn)=power_cur;
      samples.trtPre(i-numBurn)=trtPre_cur;
      samples.trtAct(i-numBurn)=trtAct_cur;
      // samples.xi(i-numBurn)=0; // don't sample here
      
      samples.ll(i-numBurn)=ll_cur;
    }
  }

  // get likelihood evaluated at posterior mean
  intcp_can = arma::mean(samples.intcp);
  for(i=0; i<covar; i++)
    beta_can(i) = arma::mean(samples.beta.col(i));
  XcovBeta_can=Xcov*beta_can;
  alpha_can = arma::mean(samples.alpha);
  alphaW_can*=arma::as_scalar(alpha_can/alpha_cur);
  power_can = arma::mean(samples.power);
  alphaW_can=alpha_can*d/arma::pow(cc,power_can);
  trtPre_can = arma::mean(samples.trtPre);
  trtAct_can = arma::mean(samples.trtAct);
  samples.llPt = ll();

  samples.Dbar=-2*arma::mean(samples.ll);
  samples.pD=samples.Dbar - -2*samples.llPt;
  samples.DIC=samples.pD + samples.Dbar;


  if(display)
    printf("\33[2K\r");
}



double GravityMcmc::ll(){
  return ll(0,1);
}



double GravityMcmc::ll(int const b, int const B){
  int i,j,k,T0,T1;
  double wontProb,prob,expProb,baseProb,baseProbInit;

  T0=(b>0 ? block2Time(b,B) : 1); // T0 must greater than 0
  T1=block2Time(b+1,B);

  for(i=T0; i<T1; i++){// loop over time interval that has changed
    llVec_can(i-1)=0;
    for(j=0; j<numNodes; j++){
      if(infHist(j,i-1)==0){// if county is susceptible get infProb
	wontProb=1.0;
	// set a base number to decrease floating point operations
	if(trtPreHist(j,i-1)==0)
	  baseProbInit=intcp_can+XcovBeta_can(j);
	else
	  baseProbInit=intcp_can+XcovBeta_can(j)+trtPre_can;

	for(k=0; k<numNodes; k++){
	  if(infHist(k,i-1)==1){// if county is infected it affects the infProb
	    // calculate infProb
	    baseProb=baseProbInit;
	    baseProb-=(j > k ? alphaW_can(j,k) : alphaW_can(k,j));

	    if(trtActHist(k,i-1)==1)
	      baseProb+=trtAct_can;

	    expProb=std::exp(-baseProb);

	    wontProb*=1-(1/(1+expProb));
	  }
	}
	
	prob=1-wontProb;

	if(!(prob>0))
	  prob=std::exp(-30.0);
	if(infHist(j,i)==0)
	  llVec_can(i-1)+=std::log(1-prob);
	else
	  llVec_can(i-1)+=std::log(prob);
      }
    }
  }
  
  return arma::sum(llVec_can);
}



int GravityMcmc::block2Time(int const b, int const B){
  int TmodB=T % B;
  return (T/B)*b + (b > TmodB ? TmodB : b);
}
