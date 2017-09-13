#ifndef LSMM_aux_hpp
#define LSMM_aux_hpp

#include <stdio.h>
#include <RcppArmadillo.h>
using namespace arma;

class LSMMfit{
public:
  LSMMfit( double alpha, vec pi1_stage1, vec pi1_stage2, vec pi1, vec b, double sigma2, double omega,
          vec omegak, vec beta, vec Lq, uword iter_times_stage1, uword iter_times_stage2,
          uword iter_times_stage3, uword iter_times_stage4){
    this -> alpha = alpha;
    this -> pi1_stage1 = pi1_stage1;
    this -> pi1_stage2 = pi1_stage2;
    this -> pi1 = pi1;
    this -> b = b;
    this -> sigma2 = sigma2;
    this -> omega = omega;
    this -> omegak = omegak;
    this -> beta = beta;
    this -> Lq = Lq;
    this -> iter_times_stage1 = iter_times_stage1;
    this -> iter_times_stage2 = iter_times_stage2;
    this -> iter_times_stage3 = iter_times_stage3;
    this -> iter_times_stage4 = iter_times_stage4;
  }

  ~LSMMfit( ){
  }

  double alpha;
  vec pi1_stage1;
  vec pi1_stage2;
  vec pi1;
  vec b;
  double sigma2;
  double omega;
  vec omegak;
  vec beta;
  vec Lq;
  uword iter_times_stage1;
  uword iter_times_stage2;
  uword iter_times_stage3;
  uword iter_times_stage4;
};


class LFMfit{
public:
  LFMfit( double alpha, vec pi1_stage1, vec pi1, vec b, uword iter_times_stage1, uword iter_times_stage2){
    this -> alpha = alpha;
    this -> pi1_stage1 = pi1_stage1;
    this -> pi1 = pi1;
    this -> b = b;
    this -> iter_times_stage1 = iter_times_stage1;
    this -> iter_times_stage2 = iter_times_stage2;
  }

  ~LFMfit( ){
  }

  double alpha;
  vec pi1_stage1;
  vec pi1;
  vec b;
  uword iter_times_stage1;
  uword iter_times_stage2;
};

class TGfit{
public:
  TGfit( double alpha, double pi1_, vec pi1, uword iter_times){
    this -> alpha = alpha;
    this -> pi1_ = pi1_;
    this -> pi1 = pi1;
    this -> iter_times = iter_times;
  }

  ~TGfit( ){
  }

  double alpha;
  double pi1_;
  vec pi1;
  uword iter_times;
};

#endif
