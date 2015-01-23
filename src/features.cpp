#include "features.hpp"


void reconstructNetwork(std::vector<int> & network,
			std::vector<int> & treated){
  int change=1,i,j,k,N=(int)treated.size();
  while(change){
    change=0;
    for(i=0; i<N; i++){
      if(treated.at(i)){
	for(j=0; j<N; j++){
	  if(network.at(i*N + j)){
	    for(k=j; k<N; k++){
	      if(network.at(i*N + k)){
		if(!network.at(j*N + k)){
		  change=1;
		  network.at(j*N + k)=network.at(k*N + j)=1;
		}
	      }
	    }
	  }
	}
      }
    }
  }
}
			

