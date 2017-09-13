#include "Main.hpp"
#include "LSMM_aux.hpp"
#include "EM_P.hpp"
#include "EM_PZ.hpp"
#include "EM_ZA.hpp"
#include "EM_PZA.hpp"

LSMMfit* Main(vec Pvalue, mat Z, mat A, uword M, double alpha, double pi1_, uword maxiter, double tol){

  uword L = Z.n_cols;
  uword K = A.n_cols;

  printf("\nWarm starts... \n \n");
  // Stage 1
  printf("Stage 1:  Two Groups Model using only p-values \n");
  vec pi1_stage1 = zeros<vec>(M);
  uword iter_times_stage1 = 0;
  EM_P(Pvalue, M, alpha, pi1_, pi1_stage1, iter_times_stage1, maxiter, tol);
  printf("  iterations=%d \n", iter_times_stage1);

  // Stage 2
  printf("Stage 2: Latent fixed-effect Model using p-values and fixed effects \n");
  vec b = zeros<vec>(L);
  b[0] = log(pi1_/(1-pi1_));
  vec pi1_stage2 = zeros<vec>(M);
  uword iter_times_stage2 = 0;
  EM_PZ(Pvalue, Z, M, L, alpha, b, pi1_stage2, iter_times_stage2, maxiter, tol);
  printf("  iterations=%d \n", iter_times_stage2);

  // Stage 3
  printf("Stage 3: logistic sparse mixed model using fixed effects, random effects and gamma(suppose gamma is known) \n");
  double sigma2 = 1;
  double omega = 0.5;
  vec omegak = zeros<vec>(K);
  vec mu = zeros<vec>(K);
  vec xi = abs(Z*b);
  uword iter_times_stage3 = 0;
  EM_ZA(Z, A, M, L, K, b, pi1_stage2, sigma2, omega, omegak, mu, xi, iter_times_stage3, maxiter, tol);
  printf("  iterations=%d \n", iter_times_stage3);


  printf("\nMain iterations... \n \n");
  // Stage 4
  printf("Main: Latent Sparse Mixed Model using p-values, fixed effects and random effects \n");
  vec pi1 = pi1_stage2;
  vec beta = zeros<vec>(K);
  vec LQ = zeros<vec>(maxiter);
  uword iter_times_stage4 = 0;
  EM_PZA(Pvalue, Z, A, M, L, K, alpha, b, pi1, sigma2, omega, omegak, mu, xi, beta, LQ, iter_times_stage4,
         maxiter, tol);
  printf("  iterations=%d \n", iter_times_stage4);
  vec Lq = LQ.subvec(0, iter_times_stage4-1);

  LSMMfit* fit = new LSMMfit(alpha, pi1_stage1, pi1_stage2, pi1, b, sigma2, omega, omegak, beta, Lq,
                           iter_times_stage1, iter_times_stage2, iter_times_stage3, iter_times_stage4);
  return fit;
};

LFMfit* Main(vec Pvalue, mat Z, uword M, double alpha, double pi1_, uword maxiter, double tol){

  uword L = Z.n_cols;

  printf("\nWarm starts... \n \n");
  // Stage 1
  printf("Stage 1:  Two Groups Model using only p-values \n");
  vec pi1_stage1 = zeros<vec>(M);
  uword iter_times_stage1 = 0;
  EM_P(Pvalue, M, alpha, pi1_, pi1_stage1, iter_times_stage1, maxiter, tol);
  printf("  iterations=%d \n", iter_times_stage1);

  printf("\nMain iterations... \n \n");
  // Stage 2
  printf("Main: Latent fixed-effect Model using p-values and fixed effects \n");
  vec b = zeros<vec>(L);
  b[0] = log(pi1_/(1-pi1_));
  vec pi1 = zeros<vec>(M);
  uword iter_times_stage2 = 0;
  EM_PZ(Pvalue, Z, M, L, alpha, b, pi1, iter_times_stage2, maxiter, tol);
  printf("  iterations=%d \n", iter_times_stage2);

  LFMfit* fit = new LFMfit(alpha, pi1_stage1, pi1, b, iter_times_stage1, iter_times_stage2);
  return fit;
};

TGfit* Main(vec Pvalue, uword M, double alpha, double pi1_, uword maxiter, double tol){

  printf("\nMain iterations... \n \n");
  printf("Main: Two Groups Model using only p-values \n");
  vec pi1 = zeros<vec>(M);
  uword iter_times = 0;
  EM_P(Pvalue, M, alpha, pi1_, pi1, iter_times, maxiter, tol);
  printf("  iterations=%d \n", iter_times);

  TGfit* fit = new TGfit(alpha, pi1_, pi1, iter_times);

  return fit;
};
