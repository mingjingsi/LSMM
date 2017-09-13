// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "Main.hpp"
#include "LSMM_aux.hpp"
#include "Rcpp_aux.hpp"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
RcppExport SEXP LSMM_Rcpp(arma::vec& Pvalue, arma::mat& Z, arma::mat& A, double alpha=0.1, double pi1_=0.1,
                         arma::uword maxiter=1e4, double tol=1e-6){
  uword M = Pvalue.n_rows;

  LSMMfit* fit = Main(Pvalue, Z, A, M, alpha, pi1_, maxiter, tol);
  return wrap_fit(fit);
}

// [[Rcpp::export]]
RcppExport SEXP LFM_Rcpp(arma::vec& Pvalue, arma::mat& Z, double alpha=0.1, double pi1_=0.1,
                          arma::uword maxiter=1e4, double tol=1e-6){
  uword M = Pvalue.n_rows;

  LFMfit* fit = Main(Pvalue, Z, M, alpha, pi1_, maxiter, tol);
  return wrap_fit_LFM(fit);
}

// [[Rcpp::export]]
RcppExport SEXP TwoGroup_Rcpp(arma::vec& Pvalue, double alpha=0.1, double pi1_=0.1,
                         arma::uword maxiter=1e4, double tol=1e-6){
  uword M = Pvalue.n_rows;

  TGfit* fit = Main(Pvalue, M, alpha, pi1_, maxiter, tol);
  return wrap_fit_TG(fit);
}
