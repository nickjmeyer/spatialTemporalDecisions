#include "m1NmOptim.hpp"


M1NmOptimTunePar::M1NmOptimTunePar(){
  mcReps = 100;
  tol = .0001;
  ss = 2.0;
}

std::vector<double> M1NmOptimTunePar::getPar() const{
  return std::vector<double> (0);
}

void M1NmOptimTunePar::putPar(const std::vector<double> & par){
}


template class M1NmOptim<System<GravityModel,GravityParam,
				GravityModel,GravityParam>,
			 RankToyAgent<ToyFeatures0<GravityModel,GravityParam>,
				      GravityModel,GravityParam>,
			 GravityModel,GravityParam>;

template <class S, class A, class M, class MP>
M1NmOptim<S,A,M,MP>::M1NmOptim(){
  name = "M1Nm";
}

template <class S, class A, class M, class MP>
void M1NmOptim<S,A,M,MP>
::optim(const S & system,
	A & agent){

  System<M,MP,M,MP> s(system.sD,system.tD,system.fD,system.dD,
		      system.modelEst,system.modelEst,
		      system.paramEst,system.paramEst);

  M1NmData<System<M,MP,M,MP>,A,M,MP> d;
  d.s = s;
  d.a = agent;
  d.r = PlainRunner<System<M,MP,M,MP>,A>();
  d.numReps = tp.mcReps;
  d.numYears = s.fD.finalT;

  std::vector<double> par = agent.tp.getPar();
  int i,dim=agent.f.numFeatures;

  gsl_vector *x, *ss;
  x = gsl_vector_alloc(dim);
  for(i=0; i<dim; i++)
    gsl_vector_set(x,i,par.at(i));
  ss=gsl_vector_alloc(dim);
  gsl_vector_set_all(ss,tp.ss);

  gsl_multimin_function minex_func;
  minex_func.n=dim;
  minex_func.f=&M1NmObj<System<M,MP,M,MP>,A,M,MP>;
  minex_func.params=&d;

  const gsl_multimin_fminimizer_type *T=
    gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *smin = NULL;
  smin=gsl_multimin_fminimizer_alloc(T,dim);
  gsl_multimin_fminimizer_set(smin,&minex_func,x,ss);

  double curSize,size=tp.tol;
  size_t iter=0;
  int status;
  
  do{
    iter++;
    status=gsl_multimin_fminimizer_iterate(smin);
    if(status)
      break;
    curSize=gsl_multimin_fminimizer_size(smin);
    status=gsl_multimin_test_size(curSize,size);

    printf("iter % d: Q() = % 16.6f  ->  [",
    	   (int)iter,smin->fval);
    for(i=0; i<(dim-1); i++)
      printf(" % 10.6f,",gsl_vector_get(smin->x,i));
    printf(" % 10.6f ]\r",gsl_vector_get(smin->x,i));
    fflush(stdout);

  }while(status == GSL_CONTINUE && iter < 100);
  // std::cout << "\033[K";

  for(i=0; i<dim; i++)
    par.at(i) = gsl_vector_get(smin->x,i);

  agent.tp.putPar(par);

  gsl_multimin_fminimizer_free(smin);
  gsl_vector_free(x);
  gsl_vector_free(ss);
}



