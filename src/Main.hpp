#ifndef Main_hpp
#define Main_hpp
#include <stdio.h>
#include <math.h>
#include "LSMM_aux.hpp"

LSMMfit* Main(vec Pvalue, mat Z, mat A, uword M, double alpha, double pi1_, uword maxiter, double tol);

LFMfit* Main(vec Pvalue, mat Z, uword M, double alpha, double pi1_, uword maxiter, double tol);

TGfit* Main(vec Pvalue, uword M, double alpha, double pi1_, uword maxiter, double tol);

#endif
