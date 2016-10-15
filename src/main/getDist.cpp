#include <iostream>
#include <cmath>

extern "C" void getDist(double * const d,
        const int * const n_,
        const int * const neigh,
        const double * const nodesX,
        const double * const nodesY,
        const int * const preAlloc){
    int i,j,k,n=(*n_);
    // std::cout << "preAlloc: " << (*preAlloc) << std::endl;
    if(!(*preAlloc)){
        for(i=0; i<n; i++){
            for(j=(i+1); j<n; j++)
                if(i!=j)
                    d[i*n + j] = d[j*n + i] = -1;
                else
                    d[i*n + j] = 0;
        }
    }
    double dIJ;
    for(i=0; i<n; i++)
        for(j=(i+1); j<n; j++)
            if(d[i*n + j] < 0 && neigh[i*n + j]){
                // dIJ =  std::sqrt(std::pow(nodesX[i] - nodesX[j],2.0) +
                // 		 std::pow(nodesY[i] - nodesY[j],2.0));
                // d[i*n + j] = d[j*n + i] = std::floor(dIJ*100000)/100000;
                d[i*n + j] = d[j*n + i] = 1.0;
            }
    int change;
    double dIK,dJK;
    do{
        change=0;
        for(i=0; i<n; i++){
            for(j=0; j<n; j++){
                if(i!=j){
                    dIJ = d[i*n + j];
                    if(dIJ >=0){
                        for(k=(i+1); k<n; k++){
                            if(j!=k){
                                dIK = d[i*n + k];
                                dJK = d[j*n + k];
                                if(dJK >=0 &&
                                        (dIK <0 || (dIJ + dJK < dIK))){
                                    d[i*n + k] = d[k*n + i] = dIJ + dJK;
                                    change=1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }while(change);
}
