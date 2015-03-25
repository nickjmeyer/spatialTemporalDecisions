#ifdef NJM_USE_MKL
#include "pardisoSymWrap.hpp"

std::vector<double> pardisoSymWrap(const std::vector<int> & iaVec,
				   const std::vector<int> & jaVec,
				   const std::vector<double> & aVec,
				   const std::vector<double> & bVec){
  /* Matrix data. */
  MKL_INT n = bVec.size();
  int i_,I_;

  I_ = iaVec.size();
  MKL_INT *ia = new MKL_INT[I_];
  for(i_ = 0; i_ < I_; ++i_)
    ia[i_] = iaVec.at(i_);
  
  I_ = jaVec.size();
  MKL_INT *ja = new MKL_INT[I_];
  for(i_ = 0; i_ < I_; ++i_)
    ja[i_] = jaVec.at(i_);

  I_ = aVec.size();
  double * a = new double[I_];
  for(i_ = 0; i_ < I_; ++i_)
    a[i_] = aVec.at(i_);

  I_ = n;
  double * b = new double[I_];
  for(i_ = 0; i_ < I_; ++i_)
    b[i_] = bVec.at(i_);

  double * x = new double[I_];
  
  MKL_INT mtype = -2;       /* Real symmetric matrix */
  /* RHS and solution vectors. */
  // double b[8], x[8];
  MKL_INT nrhs = 1;     /* Number of right hand sides. */
  /* Internal solver memory pointer pt, */
  /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
  /* or void *pt[64] should be OK on both architectures */
  void *pt[64];
  /* Pardiso control parameters. */
  MKL_INT iparm[64];
  MKL_INT maxfct, mnum, phase, error, msglvl;
  /* Auxiliary variables. */
  MKL_INT i;
  double ddum;          /* Double dummy */
  MKL_INT idum;         /* Integer dummy. */
  /* -------------------------------------------------------------------- */
  /* .. Setup Pardiso control parameters. */
  /* -------------------------------------------------------------------- */
  for ( i = 0; i < 64; i++ )
    {
      iparm[i] = 0;
    }
  iparm[0] = 1;         /* No solver default */
  iparm[1] = 2;         /* Fill-in reordering from METIS */
  iparm[3] = 0;         /* No iterative-direct algorithm */
  iparm[4] = 0;         /* No user fill-in reducing permutation */
  iparm[5] = 0;         /* Write solution into x */
  iparm[6] = 0;         /* Not in use */
  iparm[7] = 2;         /* Max numbers of iterative refinement steps */
  iparm[8] = 0;         /* Not in use */
  iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
  iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
  iparm[11] = 0;        /* Not in use */
  iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
  iparm[13] = 0;        /* Output: Number of perturbed pivots */
  iparm[14] = 0;        /* Not in use */
  iparm[15] = 0;        /* Not in use */
  iparm[16] = 0;        /* Not in use */
  iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
  iparm[18] = -1;       /* Output: Mflops for LU factorization */
  iparm[19] = 0;        /* Output: Numbers of CG Iterations */
  maxfct = 1;           /* Maximum number of numerical factorizations. */
  mnum = 1;         /* Which factorization to use. */
  msglvl = 0;           /* Print statistical information in file */
  error = 0;            /* Initialize error flag */

  /* -------------------------------------------------------------------- */
  /* .. Initialize the internal solver memory pointer. This is only */
  /* necessary for the FIRST call of the PARDISO solver. */
  /* -------------------------------------------------------------------- */
  for ( i = 0; i < 64; i++ )
    {
      pt[i] = 0;
    }
  /* -------------------------------------------------------------------- */
  /* .. Reordering and Symbolic Factorization. This step also allocates */
  /* all memory that is necessary for the factorization. */
  /* -------------------------------------------------------------------- */
  phase = 11;
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	   &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
  if ( error != 0 )
    {
      // printf ("\nERROR during symbolic factorization: %d", error);
      throw(1);
    }
  // printf ("\nReordering completed ... ");
  // printf ("\nNumber of nonzeros in factors = %d", iparm[17]);
  // printf ("\nNumber of factorization MFLOPS = %d", iparm[18]);
  /* -------------------------------------------------------------------- */
  /* .. Numerical factorization. */
  /* -------------------------------------------------------------------- */
  phase = 22;
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	   &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
  if ( error != 0 )
    {
      // printf ("\nERROR during numerical factorization: %d", error);
      throw(1);
    }
  // printf ("\nFactorization completed ... ");
  /* -------------------------------------------------------------------- */
  /* .. Back substitution and iterative refinement. */
  /* -------------------------------------------------------------------- */
  phase = 33;
  iparm[7] = 2;         /* Max numbers of iterative refinement steps. */
  // /* Set right hand side to one. */
  // for ( i = 0; i < n; i++ )
  //   {
  //     b[i] = 1;
  //   }
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	   &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
  if ( error != 0 )
    {
      // printf ("\nERROR during solution: %d", error);
      throw(1);
    }
  // printf ("\nSolve completed ... ");
  // printf ("\nThe solution of the system is: ");
  // for ( i = 0; i < n; i++ )
  //   {
  //     printf ("\n x [%d] = % f", i, x[i]);
  //   }
  // printf ("\n");

  std::vector<double> xVec;
  xVec.reserve(n);
  for(i = 0; i < n; ++i)
    xVec.push_back(x[i]);
  
  /* -------------------------------------------------------------------- */
  /* .. Termination and release of memory. */
  /* -------------------------------------------------------------------- */
  phase = -1;           /* Release internal memory. */
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	   &n, &ddum, ia, ja, &idum, &nrhs,
	   iparm, &msglvl, &ddum, &ddum, &error);

  
  delete[] ia;
  delete[] ja;
  delete[] a;
  delete[] b;
  delete[] x;
  
  return xVec;
}


void mat2Raw(const Eigen::SparseMatrix<double> & mat,
	     std::vector<int> & iaVec,
	     std::vector<int> & jaVec,
	     std::vector<double> & aVec){

  // Eigen iterates row major so rows look like columns to pardiso
  // this code has been checked....trust it
  int k,K = mat.outerSize(),n_elem = 0,col,row,p;
  p = -1;
  
  iaVec.clear();
  iaVec.reserve(mat.outerSize());
  
  jaVec.clear();
  jaVec.reserve(mat.nonZeros());
  
  aVec.clear();
  aVec.reserve(mat.nonZeros());
  
  for(k = 0; k < K; ++k){
    for(Eigen::SparseMatrix<double>::InnerIterator it(mat,k);
	it; ++it){
      col = it.col();
      row = it.row();
      if(row>=col){
	n_elem++;
	if(col > p){
	  iaVec.push_back(n_elem);
	  p = col;
	}
	jaVec.push_back(row + 1);
	aVec.push_back(it.value());
      }
    }
  }
  iaVec.push_back(n_elem + 1);
}


Eigen::VectorXd pardisoSolve(const Eigen::SparseMatrix<double> & mat,
			     const Eigen::VectorXd & vec){
  omp_set_num_threads(1);
  
  std::vector<int> ia,ja;
  std::vector<double> a,b,x;
  mat2Raw(mat,ia,ja,a);

  int i,I = vec.size();
  b.reserve(I);
  for(i = 0; i < I; ++i)
    b.push_back(vec(i));

  x = pardisoSymWrap(ia,ja,a,b);

  Eigen::VectorXd xEig(I);
  for(i = 0; i < I; ++i)
    xEig(i) = x.at(i);
  return xEig;
}
  


#endif
