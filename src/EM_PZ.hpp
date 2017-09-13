#ifndef EM_PZ_hpp
#define EM_PZ_hpp
#include <RcppArmadillo.h>
using namespace arma;

void EM_PZ(vec Pvalue, mat Z, uword M, uword L, double& alpha, vec& b, vec& pi1, uword& iter_times,
           uword maxiter, double tol);

#endif
