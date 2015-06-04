#include "nullOptim.hpp"


template class NullOptim<System<ModelTimeExpCavesGDistTrendPowCon,
				ModelTimeExpCavesGDistTrendPowCon>,
			 ProxStocGDistAgent<ModelTimeExpCavesGDistTrendPowCon>,
			 ModelTimeExpCavesGDistTrendPowCon>;


template <class S, class A, class M>
NullOptim<S,A,M>::NullOptim(){
  name = "Null";
}


template <class S, class A, class M>
void NullOptim<S,A,M>::reset(){
}


template <class S, class A, class M>
void NullOptim<S,A,M>
::optim(const S & system,
	A & agent){
}



