#ifndef Rcpp_aux_hpp
#define Rcpp_aux_hpp
#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "LSMM_aux.hpp"
using namespace arma;
using namespace Rcpp;

RcppExport SEXP wrap_fit(LSMMfit* fit){
  Rcpp::List ret;
  ret["alpha"] = fit -> alpha;
  ret["pi1.stage1"] = fit -> pi1_stage1;
  ret["pi1.stage2"] = fit -> pi1_stage2;
  ret["pi1"] = fit -> pi1;
  ret["b"] = fit -> b;
  ret["sigma2"] = fit -> sigma2;
  ret["omega"] = fit -> omega;
  ret["omegak"] = fit -> omegak;
  ret["beta"] = fit -> beta;
  ret["Lq"] = fit -> Lq;
  ret["iter_times.stage1"] = fit -> iter_times_stage1;
  ret["iter_times.stage2"] = fit -> iter_times_stage2;
  ret["iter_times.stage3"] = fit -> iter_times_stage3;
  ret["iter_times.stage4"] = fit -> iter_times_stage4;

  return ret;
}

RcppExport SEXP wrap_fit_LFM(LFMfit* fit){
  Rcpp::List ret;
  ret["alpha"] = fit -> alpha;
  ret["pi1.stage1"] = fit -> pi1_stage1;
  ret["pi1"] = fit -> pi1;
  ret["b"] = fit -> b;
  ret["iter_times.stage1"] = fit -> iter_times_stage1;
  ret["iter_times.stage2"] = fit -> iter_times_stage2;

  return ret;
}

RcppExport SEXP wrap_fit_TG(TGfit* fit){
  Rcpp::List ret;
  ret["alpha"] = fit -> alpha;
  ret["pi1_"] = fit -> pi1_;
  ret["pi1"] = fit -> pi1;
  ret["iter_times"] = fit -> iter_times;

  return ret;
}

#endif

