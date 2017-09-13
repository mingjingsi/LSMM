#ifndef EM_ZA_hpp
#define EM_ZA_hpp
#include <RcppArmadillo.h>
using namespace arma;

void EM_ZA(mat Z, mat A, uword M, uword L, uword K, vec& b, vec gamma, double& sigma2, double& omega,
           vec& omegak, vec& mu, vec& xi, uword& iter_times, uword maxiter, double tol);

#endif
