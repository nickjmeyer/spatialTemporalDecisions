#include "dataDepth.hpp"


double halfPlaneDepth(const double u, const double v, const int n,
		      const std::vector<double> & x,
		      const std::vector<double> & y){
  std::vector<double> alpha(n);
  std::vector<int> f(n);
  std::fill(alpha.begin(),alpha.end(),0);
  std::fill(f.begin(),f.end(),0);

  int numh,nt,i,j,nn,nu,ja,jb,nn2,nf,gi,ki;
  double p,p2,eps,d,xu,yu,angle,alphk,betak;
  
  numh=0;

  if(n<1)
    return 0;

  p = std::acos(-1.0);
  p2 = p*2.0;
  eps = 0.00000001;
  nt = 0;


  std::vector<double>::const_iterator xIt,yIt;
  
  for(i=0,xIt=x.begin(),yIt=y.begin(); i<n; i++,xIt++,yIt++){
    xu = (*xIt - u);
    yu = (*yIt - v);
    d = std::sqrt(xu*xu + yu*yu);

    xu /= d;
    yu /= d;
    
    if(d <= eps)
      nt++;
    else{
      if(std::abs(xu) > std::abs(yu)){
	if(*xIt >= u){
	  alpha.at(i-nt) = std::asin(yu);
	  if(alpha.at(i-nt) < 0.0)
	    alpha.at(i-nt) += p2;
	}
	else
	  alpha.at(i-nt) = p - std::asin(yu);
      }
      else{
	if(*yIt >= v)
	  alpha.at(i-nt) = std::acos(xu);
	else
	  alpha.at(i-nt) = p2 - std::acos(xu);
      }
      if(alpha.at(i-nt) >= (p2 - eps))
	alpha.at(i-nt) = 0.0;
    }
  }

  nn = n - nt;

  int pass = 0;
  if(!pass && nn <=1)
    pass = 1;

  if(!pass){
    std::vector<double>::iterator beg,end;
    beg = end = alpha.begin();
    std::advance(end,nn);
    std::sort(beg,end);

    angle = alpha.at(0) - alpha.at(nn-1) + p2;
    for(i=1; i<nn; i++)
      angle = std::max(angle, alpha.at(i) - alpha.at(i-1));
  }

  if(!pass && angle > (p + eps))
    pass = 1;

  if(!pass){
    angle = alpha.at(0);

    nu = 0;
    for(i=0; i<nn; i++){
      alpha.at(i) -= angle;
      if(alpha.at(i) < (p - eps))
	nu++;
    }
  }

  if(!pass && nu >= nn)
    pass = 1;


  if(!pass){
    ja = 1;
    jb = 1;

    alphk = alpha.at(0);
    betak = alpha.at(nu) - p;

    nn2 = nn * 2;
    i = nu;
    nf = nn;

    for(j=0; j<nn2; j++){
      if((alphk + eps) < betak){
	nf++;

	if(ja < nn){
	  alphk = alpha.at(ja++);
	}
	else
	  alphk = p2 + 1.0;
      }
      else{
	i++;
	if(i == (nn+1)){
	  i = 1;
	  nf -= nn;
	}
	f.at(i-1) = nf;

	if(jb < nn){
	  jb++;
	  if((jb + nu) <= nn)
	    betak = alpha.at(jb + nu - 1) - p;
	  else
	    betak = alpha.at(jb + nu - nn - 1) + p;
	}
	else
	  betak = p2 + 1.0;
      }
    }

    gi = 0;
    ja = 1;

    angle = alpha.at(0);

    numh = std::min(f.at(0), nn - f.at(0));

    for(i=1; i<nn; i++){
      if(alpha.at(i) <= (angle + eps))
	ja++;
      else{
	gi += ja;
	ja = 1;
	angle = alpha.at(i);
      }
      ki = f.at(i) - gi;
      numh = std::min(numh, std::min(ki, nn - ki));
    }
  }

  numh += nt;
  return ((double)numh) / ((double) n);
}


