#include "densityEst.hpp"

DensityEst::DensityEst(const std::vector<double> xi)
  : xi(xi), h(optH(xi)){
}




double DensityEst::gKernel(const double x){
  return std::exp(-x*x/2)/std::sqrt(2*M_PI);
}

double DensityEst::gKernelP(const double x){
  return -x*std::exp(-x*x/2)/std::sqrt(2*M_PI);
}

double DensityEst::optH(const std::vector<double> x){
  const int n = int(x.size());
  const double mn = std::accumulate(x.begin(),x.end(),0.0)/double(n);
  const double var = std::accumulate(x.begin(),x.end(),0.0,
				     [mn,n](const double a, const double b){
				       return a + (b-mn)*(b-mn);
				     })/double(n-1);
  const double sd = std::sqrt(var);

  return std::pow(4.0/(double(n)*3.0),0.2)*sd;
}

double DensityEst::eval(const double x){
  return std::accumulate(xi.begin(),xi.end(),0.0,
			 [x,this](const double a, const double b){
			   return a + gKernel((x-b)/h);
			 })/(double(xi.size())*h);
}

double DensityEst::deriv(const double x){
  return std::accumulate(xi.begin(),xi.end(),0.0,
			 [x,this](const double a, const double b){
			   return a + gKernelP((x-b)/h);
			 })/(double(xi.size())*h*h);
}


double DensityEst::f(const gsl_vector * x, void * params){
  DensityEst * de = static_cast<DensityEst*>(params);
  return -de->eval(gsl_vector_get(x,0));
}

void DensityEst::df(const gsl_vector * x, void * params,
		    gsl_vector * g){
  DensityEst * de = static_cast<DensityEst*>(params);
  gsl_vector_set(g,0,-de->deriv(gsl_vector_get(x,0)));
}

void DensityEst::fdf(const gsl_vector * x, void * params,
		     double * f, gsl_vector * g){
  *f = DensityEst::f(x,params);
  DensityEst::df(x,params,g);
}




std::pair<double,double> DensityEst::max(){
  std::vector<double> xiSorted = xi;
  std::sort(xiSorted.begin(),xiSorted.end());

  std::priority_queue<std::pair<double,double> > ps;

  const int nLen = int(xi.size());
  const int nTry = std::log(double(nLen)) + 1.0;
  int n;
  for(n = 0; n < nTry; ++n){
    size_t iter = 0;
    int status;
    gsl_vector *x;
    int dim = 1;
    x = gsl_vector_alloc(dim);
    gsl_vector_set(x,0,xiSorted.at(n*nLen/nTry));
    gsl_multimin_function_fdf minex_func;

    minex_func.n = dim;
    minex_func.f = &DensityEst::f;
    minex_func.df = &DensityEst::df;
    minex_func.fdf = &DensityEst::fdf;
    minex_func.params = this;

    const gsl_multimin_fdfminimizer_type * T;
    gsl_multimin_fdfminimizer * s;
    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    s = gsl_multimin_fdfminimizer_alloc(T,dim);

    gsl_multimin_fdfminimizer_set(s,&minex_func, x, 0.01, 1e-6);


    do {
      ++iter;
      status = gsl_multimin_fdfminimizer_iterate(s);

      if(status){
	break;
      }

      status = gsl_multimin_test_gradient(s->gradient,1e-6);


    } while(status == GSL_CONTINUE && iter < 100);

    std::pair<double,double> res(-s->f,gsl_vector_get(s->x,0));

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);

    ps.push(res);
  }

  return ps.top();
}
