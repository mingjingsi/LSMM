#ifndef EM_PZA_hpp
#define EM_PZA_hpp
#include <RcppArmadillo.h>
using namespace arma;

void EM_PZA(vec Pvalue, mat Z, mat A, uword M, uword L, uword K, double& alpha, vec& b, vec& pi1,
            double& sigma2, double& omega, vec& omegak, vec& mu, vec& xi, vec&beta, vec& LQ,
            uword& iter_times, uword maxiter, double tol);

#endif

