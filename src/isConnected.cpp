#include <cmath>

extern "C" void isConnected(int * const b,
			    const int * const n_,
			    const int * const neigh){
  int i,j,t=0,tOld,n=(*n_),c[n],d[n];
  for(i=0; i<n; i++){
    d[i]=0;
    if(neigh[i]){
      c[i]=1;
      t++;
    }
    else
      c[i]=0;
  }
  do{
    tOld=t;
    for(i=0; i<n; i++)
      if(c[i] && !d[i]){
	d[i]=1;
	for(j=0; j<n; j++)
	  if(neigh[i*n + j] && !c[j]){
	    c[j]=1;
	    t++;
	  }
      }
  }while(tOld!=t);
  *b = (t == n);
}


  
			    
