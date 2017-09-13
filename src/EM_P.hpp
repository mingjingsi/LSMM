#ifndef EM_P_hpp
#define EM_P_hpp
#include <RcppArmadillo.h>
using namespace arma;

void EM_P(vec Pvalue, uword M, double& alpha, double& pi1_, vec& pi1, uword& iter_times,
          uword maxiter, double tol);

#endif
