#include "sample.hpp"

template<typename T>
void sampleVals(std::vector<T> * const sample,
		std::vector<T> const * const values,
		std::vector<double> * const probs,
		int const sampleSize){
  sample->clear();
  sample->reserve(sampleSize);

  std::vector<double> probsCpy = *probs;
  int i,j,k,numValues=(int)values->size();
  double r, probSum;

  for(i=0; i<sampleSize; i++){
    r=njm::runif01();
    j=0;
    probSum=0;
    while(probSum<=r && j < numValues){ // 0 <= r < 1
      probSum+=probsCpy.at(j);
      j++;
    }
    j--;
    sample->push_back(values->at(j));

    for(k=0; k<numValues; k++)
      if(k!=j)
	probsCpy.at(k)/=(1.0-probsCpy.at(j));

    probsCpy.at(j)=0.0;

  }
}

// templates
template void sampleVals(std::vector<int> * const sample,
			 std::vector<int> const * const values,
			 std::vector<double> * const probs,
			 int const sampleSize);

template void sampleVals(std::vector<double> * const sample,
			 std::vector<double> const * const values,
			 std::vector<double> * const probs,
			 int const sampleSize);
////////////////////////////////////////////////////////////////////////////////




template<typename T>
void sampleVals(std::vector<T> * const sample,
		std::vector<T> const * const values,
		int const sampleSize){
  sample->clear();
  sample->reserve(sampleSize);
  
  std::vector<double> probs;
  int i;
  double numValues =(double)values->size();

  for(i=0; i<numValues; i++)
    probs.push_back(1.0/numValues);

  sampleVals(sample,values,&probs,sampleSize);
}

// templates
template void sampleVals(std::vector<int> * const sample,
			 std::vector<int> const * const values,
			 int const sampleSize);

template void sampleVals(std::vector<double> * const sample,
			 std::vector<double> const * const values,
			 int const sampleSize);
////////////////////////////////////////////////////////////////////////////////


template<typename T>
void shuffle(std::vector<T> * const vec){
  std::vector<int> cpy=*vec;
  int size=vec->size();
  sampleVals(vec,&cpy,size);
}

// templates
template void shuffle(std::vector<int> * const vec);
