#include "incremAgent.hpp"


template class
IncremAgent<ModelTimeExpCavesGDistTrendPowCon,
	    ProximalGDistAgent<ModelTimeExpCavesGDistTrendPowCon>,
	    NullOptim<System<ModelTimeExpCavesGDistTrendPowCon,
			     ModelTimeExpCavesGDistTrendPowCon>,
		      ProximalGDistAgent<ModelTimeExpCavesGDistTrendPowCon>,
		      ModelTimeExpCavesGDistTrendPowCon> >;



template <class M, class A, class O>
IncremAgent<M,A,O>::IncremAgent(){
  tp.N = 100;
  tp.mcReps = 10;
  
  name="increm";
}


template <class M, class A, class O>
void IncremAgent<M,A,O>::reset(){
}

  
template <class M, class A, class O>
void IncremAgent<M,A,O>::applyTrt(const SimData & sD,
				  TrtData & tD,
				  const FixedData & fD,
				  const DynamicData & dD,
				  M & m){
  if(sD.notInfec.empty())
    return;

  System<M,M> s(sD,tD,fD,dD,m,m);

  O o;
  A a;
  o.optim(s,a);

  numAct = getNumAct(sD,tD,fD,dD);
  numPre = getNumPre(sD,tD,fD,dD);

  int numTot = numAct + numPre;

  const int preVal = 0, actVal = 1;

  std::priority_queue<std::pair<double,int> > pq;

  int i;
  for(i = 0; i < numPre; ++i)
    pq.push(std::pair<double,int>(njm::runif01(),preVal));
  for(i = 0; i < numAct; ++i)
    pq.push(std::pair<double,int>(njm::runif01(),actVal));

  std::vector<int> infCand = sD.infected;
  std::vector<int> notCand = sD.notInfec;
  std::vector<int>::const_iterator it,beg,end;
  int type;
  for(i = 0; i < numTot; ++i){
    printf("Round: %3d\n",i);
    
    type = pq.top().second;
    pq.pop();

    // setup iterator bounds
    if(type == preVal){
      beg = notCand.begin();
      end = notCand.end();
    }
    else if(type == actVal){
      beg = infCand.begin();
      end = infCand.end();
    }
    else{
      std::cout << "In IncremAgent::applyTrt : invalid type"
		<< std::endl;
      throw(1);
    }

    // eval candidates
    typedef std::pair<double,std::vector<int>::const_iterator> DblIt;
    std::priority_queue<DblIt> res;
    for(it = beg; it != end; ++it){
      printf("Test: %3d\n",*it);
      
      if(type == preVal)
	s.tD_r.p.at(*it) = 1;
      else if(type == actVal)
	s.tD_r.a.at(*it) = 1;
      
      res.push(DblIt(-eval(s,a),it));
      
      if(type == preVal)
	s.tD_r.p.at(*it) = 0;
      else if(type == actVal)
	s.tD_r.a.at(*it) = 0;
    }

    // treat the best and remove from candidates
    std::vector<int>::const_iterator best = res.top().second;
    if(type == preVal){
      tD.p.at(*best) = 1;
      s.tD_r.p.at(*best) = 1;
      notCand.erase(best);
    }
    else if(type == actVal){
      tD.a.at(*best) = 1;
      s.tD_r.a.at(*best) = 1;
      infCand.erase(best);
    }
  }


#ifndef NJM_NO_DEBUG
  int totPre = 0,totAct = 0;
  // check if valid treatments are given to valid locations
  for(i = 0; i < fD.numNodes; i++){
    if(tD.p.at(i) != 1 && tD.p.at(i) != 0){
      std::cout << "Prevenative treatment not 1 or 0"
		<< ": " << tD.p.at(i)
		<< std::endl;
      throw(1);
    }
    else if(tD.a.at(i) != 1 && tD.a.at(i) != 0){
      std::cout << "Active treatment not 1 or 0"
		<< std::endl;
      throw(1);
    }
    else if(tD.a.at(i) == 1 && sD.status.at(i) < 2){
      std::cout << "Not infected receiving active treatment"
		<< std::endl;
      throw(1);
    }
    else if(tD.p.at(i) == 1 && sD.status.at(i) >= 2){
      std::cout << "Infected receiving preventative treament"
		<< std::endl;
      throw(1);
    }
    else if(tD.a.at(i) == 1)
      totAct++;
    else if(tD.p.at(i) == 1)
      totPre++;
  }

  // check if total number of treatments are correct
  if(totAct != numAct){
    std::cout << "Not correct amount of active treatments."
	      << std::endl
	      << "Should be " << numAct << " but is " << totAct << "."
	      << std::endl
	      << "Number of infected nodes is " << sD.numInfected
	      << "(" << sD.infected.size() << ")"
	      << std::endl;
    throw(1);
  }
  else if(totPre != numPre){
    std::cout << "Not correct amount of preventative treatments."
	      << std::endl
	      << "Should be " << numPre << " but is " << totPre << "."
	      << std::endl
	      << "Number of not infected nodes is " << sD.numNotInfec
	      << "(" << sD.notInfec.size() << ")"
	      << std::endl;
    throw(1);
  }
#endif
}


template <class M, class A, class O>
double IncremAgent<M,A,O>::eval(System<M,M> s,
				A a){
  PlainRunner<System<M,M>,A> r;
  printf("eval: %3d\n",tp.N);

  int i;
  double tot = 0.0;
  for(i = 0; i < tp.N; ++i){
    s.revert();

    int sumPre, sumAct;
    sumPre = std::accumulate(s.tD.p.begin(),s.tD.p.end(),0);
    sumAct = std::accumulate(s.tD.a.begin(),s.tD.a.end(),0);

    printf("(%3d,%3d)\n",sumPre,sumAct);

    s.nextPoint();
    
    tot += r.run(s,a,tp.mcReps,s.fD.finalT).smean();
  }
  return tot/(double(tp.N));
}



std::vector<double> IncremTuneParam::getPar() const {
  return std::vector<double>();
}



void IncremTuneParam::putPar(const std::vector<double> & par){
}



