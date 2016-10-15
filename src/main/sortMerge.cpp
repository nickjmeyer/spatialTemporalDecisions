#include "sortMerge.hpp"


template<typename T>
void sortByValue(std::vector<T> * A){
    int n=A->size();
    int width,i;
    std::vector<T> B;
    for(i=0; i<n; i++)
        B.push_back(0);
    for(width=1; width<n; width*=2){ // increase run-width by factor of 2
        for(i=0; i<n; i+=(2*width)) // sort and merge run
            mergeByValue(A,&B,i,std::min(i+width,n), std::min(i+2*width,n));
        *A=B;
    }
}

template void sortByValue(std::vector<int> * const A);
template void sortByValue(std::vector<double> * const A);



template<typename T>
void mergeByValue(std::vector<T> const * const A,
        std::vector<T> * const B,
        int const iLeft, int const iRight, int const iEnd){
    int i0=iLeft;
    int i1=iRight;
    int j;
    for(j=iLeft; j<iEnd; j++){
        if(i0<iRight &&
                (i1 >= iEnd || A->at(i0) < A->at(i1))
                ){ // check bounds and compare elements
            B->at(j)=A->at(i0);
            i0++;
        }
        else{
            B->at(j)=A->at(i1);
            i1++;
        }
    }
}

template void mergeByValue(std::vector<int> const * const A,
        std::vector<int> * const B,
        int const iLeft, int const iRight, int const iEnd);
template void mergeByValue(std::vector<double> const * const A,
        std::vector<double> * const B,
        int const iLeft, int const iRight, int const iEnd);



void sortByVec(std::vector<int> * A,
        std::vector<double> * const int2double){
    int n=A->size();
    if(n!=(int)int2double->size()){
        std::cout << "In sortByVec: invalid lengths\n";
        throw(1);
    }

    std::map<int,double> int2doubleMap;
    int i;
    for(i=0; i<n; i++)
        int2doubleMap[A->at(i)]=int2double->at(i);

    sortByMap(A,&int2doubleMap);
    sortByValue(int2double);
}




void sortByMap(std::vector<int> * A,
        std::map<int,double> const * const int2double){
    int n=A->size();
    int width,i;
    std::vector<int> B;
    B.reserve(n);
    for(i=0; i<n; i++)
        B.push_back(0);
    for(width=1; width<n; width*=2){ // increase run-width by factor of 2
        for(i=0; i<n; i+=(2*width)) // sort and merge run
            mergeByMap(A,int2double,&B,
                    i,std::min(i+width,n), std::min(i+2*width,n));
        *A=B;
    }
}



void mergeByMap(std::vector<int> const * const A,
        std::map<int,double> const * const int2double,
        std::vector<int> * const B,
        int const iLeft, int const iRight, int const iEnd){
    int i0=iLeft;
    int i1=iRight;
    int j;
    for(j=iLeft; j<iEnd; j++){
        if(i0<iRight &&
                (i1 >= iEnd || (int2double->at(A->at(i0)) < int2double->at(A->at(i1))))
                ){ // check bounds and compare elements
            B->at(j)=A->at(i0);
            i0++;
        }
        else{
            B->at(j)=A->at(i1);
            i1++;
        }
    }
}
