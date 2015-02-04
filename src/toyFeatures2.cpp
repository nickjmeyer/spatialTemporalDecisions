#include "toyFeatures2.hpp"

std::vector<double> ToyFeatures2TuneParam::getPar() const{
  return std::vector<double>(0);
}


void ToyFeatures2TuneParam::putPar(const std::vector<double> & par){
  // do nothing!
}





template class ToyFeatures2<GravityModel,GravityParam>;

template class ToyFeatures2<RangeModel,RangeParam>;

template class ToyFeatures2<EbolaModel,EbolaParam>;




template <class M, class MP>
int ToyFeatures2<M,MP>::numFeatures = 4;


template <class M, class MP>
void ToyFeatures2<M,MP>::preCompData(const SimData & sD,
				     const TrtData & tD,
				     const FixedData & fD,
				     const DynamicData & dD,
				     const M & m,
				     MP & mP){
  // pre compute stuff

  // ////////////////////////////////////////
  // // temporary timing stuff
  // static int numTimers=3;
  // static std::vector<std::chrono::milliseconds> preTime(numTimers);
  // static int is_init=0;
  // if(!is_init)
  //   for(int i=0; i<numTimers; i++)
  //     preTime.at(i) = std::chrono::milliseconds::zero();
  // is_init=1;

  // std::chrono::time_point< std::chrono::high_resolution_clock > tick,tock;
  // std::chrono::milliseconds diff;
  // static int timeReps = 0;
  // timeReps++;

  // int curTimer=0;

  // ////////////////////////////////////////

  

  // ////////////////////////////////////////
  // tick=std::chrono::high_resolution_clock::now();
  // ////////////////////////////////////////
    
  m.load(sD,tD,fD,dD,mP);

  // ////////////////////////////////////////
  // tock=std::chrono::high_resolution_clock::now();
  // diff=std::chrono::milliseconds::zero();
  // diff += std::chrono::duration_cast< std::chrono::milliseconds >
  //   (tock.time_since_epoch());
  // diff -= std::chrono::duration_cast< std::chrono::milliseconds >
  //   (tick.time_since_epoch());
  // preTime.at(curTimer++)+=diff;
  // ////////////////////////////////////////


  
  // ////////////////////////////////////////
  // tick=std::chrono::high_resolution_clock::now();
  // ////////////////////////////////////////
  
  int i;
  subGraphNotInfec.resize(sD.numNotInfec);
  for(i=0; i<sD.numNotInfec; i++)
    subGraphNotInfec(i)=fD.subGraph.at(sD.notInfec.at(i));

  // ////////////////////////////////////////
  // tock=std::chrono::high_resolution_clock::now();
  // diff=std::chrono::milliseconds::zero();
  // diff += std::chrono::duration_cast< std::chrono::milliseconds >
  //   (tock.time_since_epoch());
  // diff -= std::chrono::duration_cast< std::chrono::milliseconds >
  //   (tick.time_since_epoch());
  // preTime.at(curTimer++)+=diff;
  // ////////////////////////////////////////
  

  // ////////////////////////////////////////
  // tick=std::chrono::high_resolution_clock::now();
  // ////////////////////////////////////////

  
  notNeigh.resize(sD.numNotInfec);
  notNeighOf.resize(sD.numNotInfec);
  std::fill(notNeigh.begin(),notNeigh.end(),
	    std::vector<std::pair<int,double> >(0));
  std::fill(notNeighOf.begin(),notNeighOf.end(),
	    std::vector<std::pair<int,int> >(0));

  notNeighNum.resize(sD.numNotInfec);
  notNeighOfNum.resize(sD.numNotInfec);
  std::fill(notNeighNum.begin(),notNeighNum.end(),0);
  std::fill(notNeighOfNum.begin(),notNeighOfNum.end(),0);
  
  
  std::vector<int>::const_iterator itD0,itD1,beg;
  int j;
  beg=sD.notInfec.begin();
  for(i=0,itD0=beg; i<sD.numNotInfec; i++,itD0++){
    for(j=0,itD1=beg; j<sD.numNotInfec; j++,itD1++)
      if(i!=j && fD.network.at((*itD0)*fD.numNodes + (*itD1))){
	notNeigh.at(i).push_back(std::pair<int,double>
				 (j,m.oneOnOne(*itD0,*itD1,sD,tD,fD,dD,mP)));
	notNeighOf.at(j).push_back(std::pair<int,int>(i,notNeighNum.at(i)));
				   
	notNeighNum.at(i)++;
	notNeighOfNum.at(j)++;
      }
  }

  // ////////////////////////////////////////
  // tock=std::chrono::high_resolution_clock::now();
  // diff=std::chrono::milliseconds::zero();
  // diff += std::chrono::duration_cast< std::chrono::milliseconds >
  //   (tock.time_since_epoch());
  // diff -= std::chrono::duration_cast< std::chrono::milliseconds >
  //   (tick.time_since_epoch());
  // preTime.at(curTimer++)+=diff;
  // ////////////////////////////////////////

  
  // ////////////////////////////////////////
  // double total;
  // printf(" preTime: ");
  // for(i=0; i<numTimers; i++){
  //   total = (double(preTime.at(i).count()))/(double(timeReps));
  //   printf("% 20.10f ",total);
  // }
  // printf("\n");
  // fflush(stdout);
  // std::cout << std::string(40,'-') << std::endl	<< std::endl;
  // ////////////////////////////////////////
  
}



template <class M, class MP>
void ToyFeatures2<M,MP>::getFeatures(const SimData & sD,
				     const TrtData & tD,
				     const FixedData & fD,
				     const DynamicData & dD,
				     const M & m,
				     MP & mP){
  // ////////////////////////////////////////
  // // temporary timing stuff
  // static std::vector<std::chrono::milliseconds> featTime(numFeatures);
  // static int is_init=0;
  // if(!is_init)
  //   for(int i=0; i<numFeatures; i++)
  //     featTime.at(i) = std::chrono::milliseconds::zero();
  // is_init=1;

  // std::chrono::time_point< std::chrono::high_resolution_clock > tick,tock;
  // std::chrono::milliseconds diff;
  // static int timeReps = 0;
  // timeReps++;
  
  // ////////////////////////////////////////
  
  
  // clear containers
  infFeat.zeros(sD.numInfected,numFeatures);
  notFeat.zeros(sD.numNotInfec,numFeatures);


  
  // start feature construction


  
  int i,j,featNum=0;
  std::vector<int>::const_iterator itD0,itD1,beg;

  // ////////////////////////////////////////
  // tick=std::chrono::high_resolution_clock::now();
  // ////////////////////////////////////////
  
  // feature 0
  infFeat.col(featNum) = 1 - arma::prod(mP.infProbsSep,1);
  notFeat.col(featNum) = 1 - arma::prod(mP.infProbsSep,0).t();
  
  // ////////////////////////////////////////
  // tock=std::chrono::high_resolution_clock::now();
  // diff=std::chrono::milliseconds::zero();
  // diff += std::chrono::duration_cast< std::chrono::milliseconds >
  //   (tock.time_since_epoch());
  // diff -= std::chrono::duration_cast< std::chrono::milliseconds >
  //   (tick.time_since_epoch());
  // featTime.at(featNum)+=diff;
  // ////////////////////////////////////////

  
  featNum++;


  // ////////////////////////////////////////
  // tick=std::chrono::high_resolution_clock::now();
  // ////////////////////////////////////////

  
  
  // feature 1
  int num;
  double modProbTot,modProb;
  std::pair<int,double> neigh;
  for(i=0; i<sD.numNotInfec; i++){
    modProbTot=0;
    num = notNeighNum.at(i);
    for(j=0; j<num; j++){
      neigh=notNeigh.at(i).at(j);
      
      modProb = 1.0 - notFeat(neigh.first,0);
      modProb *= 1.0/(1.0 + std::exp(neigh.second));
      modProbTot += 1.0 - modProb;
    }
    notFeat(i,featNum) = modProbTot*notFeat(i,0);
  }
  
  infFeat.col(featNum) = (1.0-mP.infProbsSep) * notFeat.col(featNum);


  // ////////////////////////////////////////
  // tock=std::chrono::high_resolution_clock::now();
  // diff=std::chrono::milliseconds::zero();
  // diff += std::chrono::duration_cast< std::chrono::milliseconds >
  //   (tock.time_since_epoch());
  // diff -= std::chrono::duration_cast< std::chrono::milliseconds >
  //   (tick.time_since_epoch());
  // featTime.at(featNum)+=diff;
  // ////////////////////////////////////////

  
  
  featNum++;


  // ////////////////////////////////////////
  // tick=std::chrono::high_resolution_clock::now();
  // ////////////////////////////////////////

  
  
  // feature 2
  notFeat.col(featNum) = notFeat.col(0) % subGraphNotInfec;

  infFeat.col(featNum) = (1.0 - mP.infProbsSep) * notFeat.col(0);


  // ////////////////////////////////////////
  // tock=std::chrono::high_resolution_clock::now();
  // diff=std::chrono::milliseconds::zero();
  // diff += std::chrono::duration_cast< std::chrono::milliseconds >
  //   (tock.time_since_epoch());
  // diff -= std::chrono::duration_cast< std::chrono::milliseconds >
  //   (tick.time_since_epoch());
  // featTime.at(featNum)+=diff;
  // ////////////////////////////////////////



    
    
  featNum++;

  
  // ////////////////////////////////////////
  // tick=std::chrono::high_resolution_clock::now();
  // ////////////////////////////////////////
  

  // feature 3
  std::vector<int>::const_iterator itD2,itD3;
  std::priority_queue<double> p; 
  itD2 = sD.notInfec.begin();
  itD3 = sD.infected.begin();
  double totalDist;
  for(i=0,itD0 = itD2; i<sD.numNotInfec; i++,itD0++){
    totalDist=0;
    for(j=0,itD1 = itD3; j<sD.numInfected; j++,itD1++)
      totalDist += fD.expInvDistSD.at((*itD0)*fD.numNodes + *itD1);
    totalDist /= fD.numNodes*fD.numNodes*fD.invDistSD;
    notFeat(i,featNum) = std::log(1.0+totalDist);
  }

  for(i=0,itD0 = itD3; i<sD.numInfected; i++,itD0++){
    totalDist=0;
    for(j=0,itD1 = itD2; j<sD.numNotInfec; j++,itD1++)
      totalDist += fD.expInvDistSD.at((*itD0)*fD.numNodes + *itD1);
    totalDist /= fD.numNodes*fD.numNodes*fD.invDistSD;
    infFeat(i,featNum) = std::log(1.0+totalDist);
  }



  // ////////////////////////////////////////
  // tock=std::chrono::high_resolution_clock::now();
  // diff=std::chrono::milliseconds::zero();
  // diff += std::chrono::duration_cast< std::chrono::milliseconds >
  //   (tock.time_since_epoch());
  // diff -= std::chrono::duration_cast< std::chrono::milliseconds >
  //   (tick.time_since_epoch());
  // featTime.at(featNum)+=diff;
  // ////////////////////////////////////////


  
  // ////////////////////////////////////////
  // double total;
  // std::cout << "  infMean: " << arma::mean(infFeat)
  // 	    << "    infSD: " << arma::stddev(infFeat)
  // 	    << "  notMean: " << arma::mean(notFeat)
  // 	    << "    notSD: " << arma::stddev(notFeat)
  // 	    << "     time: ";
  // for(i=0; i<numFeatures; i++){
  //   total = (double(featTime.at(i).count()))/(double(timeReps));
  //   printf("% 20.10f ",total);
  // }
  // printf("\n");
  // fflush(stdout);
  // std::cout << std::string(40,'-') << std::endl	<< std::endl;
  // ////////////////////////////////////////

  
  featNum++;

  tDPre = tD;

  if(featNum != numFeatures){
    std::cout << "Error: in getFeatures: featNum != numFeatures"
	      << std::endl;
    throw 1;
  }
    
}



template <class M, class MP>
void ToyFeatures2<M,MP>::updateFeatures(const SimData & sD,
					const TrtData & tD,
					const FixedData & fD,
					const DynamicData & dD,
					const M & m,
					MP & mP){

  // ////////////////////////////////////////
  // // temporary timing stuff
  // static std::chrono::milliseconds updTime=std::chrono::milliseconds::zero();

  // std::chrono::time_point< std::chrono::high_resolution_clock > tick,tock;
  // std::chrono::milliseconds diff;
  // static int timeReps = 0;
  // timeReps++;
  
  // ////////////////////////////////////////


  
  // ////////////////////////////////////////
  // tick=std::chrono::high_resolution_clock::now();
  // ////////////////////////////////////////
  
  int i,j,node0,now,pre,num;
  std::pair<int,int> neighOf;

  // update not infected probabilities
  for(i=0; i<sD.numNotInfec; i++){
    node0 = sD.notInfec.at(i);
    now = tD.p.at(node0);
    pre = tDPre.p.at(node0);
    
    if(now != pre && pre == 1){ // adding trt
      mP.infProbsBase.col(i) -= mP.trtPre;
      mP.setCol(i);

      num=notNeighOfNum.at(i);
      for(j=0; j<num; j++){
	neighOf = notNeighOf.at(i).at(j);

	notNeigh.at(neighOf.first).at(neighOf.second).second -= mP.trtPre;
      }
    }
    
    else if(now != pre && pre == 0){ // removing trt
      mP.infProbsBase.col(i) += mP.trtPre;
      mP.setCol(i);

      num=notNeighOfNum.at(i);
      for(j=0; j<num; j++){
	neighOf = notNeighOf.at(i).at(j);

	notNeigh.at(neighOf.first).at(neighOf.second).second += mP.trtPre;
      }
    }
  }

  // update infected probabilities
  for(i=0; i<sD.numInfected; i++){
    node0 = sD.infected.at(i);
    now = tD.a.at(node0);
    pre = tDPre.a.at(node0);
    
    if(now != pre && pre == 1){ // adding trt
      mP.infProbsBase.row(i) -= mP.trtAct;
      mP.setRow(i);
    }
    else if(now != pre && pre == 0){ // removing trt
      mP.infProbsBase.row(i) += mP.trtPre;
      mP.setRow(i);
    }
  }


  // ////////////////////////////////////////
  // tock=std::chrono::high_resolution_clock::now();
  // diff=std::chrono::milliseconds::zero();
  // diff += std::chrono::duration_cast< std::chrono::milliseconds >
  //   (tock.time_since_epoch());
  // diff -= std::chrono::duration_cast< std::chrono::milliseconds >
  //   (tick.time_since_epoch());
  // updTime+=diff;
  // ////////////////////////////////////////

  
  // ////////////////////////////////////////
  // double total;
  // total = (double(updTime.count()))/(double(timeReps));
  // printf("  updTime: % 20.10f \n",total);
  // fflush(stdout);
  // std::cout << std::string(40,'-') << std::endl	<< std::endl;
  // ////////////////////////////////////////


  getFeatures(sD,tD,fD,dD,m,mP);
}


